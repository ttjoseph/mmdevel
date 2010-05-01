#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include <mpi.h>

#define COULOMB 332.0636
#define EEL_14_SCALING_RECIP (1 / 1.2)
#define VDW_14_SCALING_RECIP (1 / 2.0)
// Flags for BondType
#define BOND 1
#define ANGLE 2
#define DIHEDRAL 4

// Loads stuff from disk
// Assumes the array is in host endianness!
int* loadIntArray(FILE *fp, int *size) {
  // Get number of elements in array
  fread(size, sizeof(int), 1, fp);
  int *data = malloc(sizeof(int) * *size);
  fread(data, sizeof(int), *size, fp);
  return data;
}

uint8_t* loadByteArray(FILE *fp, int *size) {
  // Get number of elements in array
  fread(size, sizeof(int), 1, fp);
  uint8_t *data = malloc(sizeof(uint8_t) * *size);
  fread(data, sizeof(uint8_t), *size, fp);
  return data;
}

float* loadFloatArray(FILE *fp, int *size) {
  // Get number of elements in array
  fread(size, sizeof(int), 1, fp);
  float *data = malloc(sizeof(float) * *size);
  fread(data, sizeof(float), *size, fp);
  return data;
}

int Rank, NumNodes; // For MPI

// Synchronizes data by broadcasting it from the rank 0 node
void broadcastIntArray(int **buf, int *size) {
  MPI_Bcast(size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  if(Rank > 0) // If we're not the master node, we need to allocate this buffer
    *buf = malloc(sizeof(int) * *size);
  MPI_Bcast(*buf, *size, MPI_INTEGER, 0, MPI_COMM_WORLD);
}

void syncByteArray(uint8_t **buf, int *size) {
  MPI_Bcast(size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  if(Rank > 0) // If we're not the master node, we need to allocate this buffer
    *buf = malloc(sizeof(uint8_t) * *size);
  MPI_Bcast(*buf, *size, MPI_BYTE, 0, MPI_COMM_WORLD);
}

void syncFloatArray(float **buf, int *size) {
  MPI_Bcast(size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  if(Rank > 0) // If we're not the master node, we need to allocate this buffer
    *buf = malloc(sizeof(float) * *size);
  MPI_Bcast(*buf, *size, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void bomb(char *msg) {
  printf("[%d] ERROR, aborting! %s\n", Rank, msg);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}

int Ntypes, Natoms, Nresidues;
int *NBIndices, NumNBIndices;
int *AtomTypeIndices, NumAtomTypeIndices;
float *LJ12; int NumLJ12;
float *LJ6; int NumLJ6;
float *Charges; int NumCharges;
uint8_t *BondType; int NumBondType;
int *ResidueMap, NumResidueMap;

// Synchronize parameters for the molecule
void syncMolecule() {
  MPI_Bcast(&Natoms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nresidues, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ntypes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  broadcastIntArray(&NBIndices, &NumNBIndices);
  broadcastIntArray(&AtomTypeIndices, &NumAtomTypeIndices);
  broadcastIntArray(&ResidueMap, &NumResidueMap);
  syncFloatArray(&LJ12, &NumLJ12);
  syncFloatArray(&LJ6, &NumLJ6);
  syncFloatArray(&Charges, &NumCharges);
  syncByteArray(&BondType, &NumBondType);
}

// van der Waals calculation by Lennard-Jones 12-6 method as used by AMBER. No cutoffs etc
double LennardJones(float *coords, double *decomp) {
  double energy = 0.0;
  int atom_i;
	
  for(atom_i = 0; atom_i < Natoms; atom_i++) {
    int offs_i = atom_i * 3;
    float x0 = coords[offs_i], y0 = coords[offs_i+1], z0 = coords[offs_i+2];
    // Pulled some of the matrix indexing out of the inner loop
    int nbparm_offs_i = Ntypes * (AtomTypeIndices[atom_i] - 1);
    int bondtype_offs_i = atom_i * Natoms;
    int i_res = ResidueMap[atom_i]; // Residue of atom i

    int atom_j;
    for(atom_j = 0; atom_j < atom_i; atom_j++) {
      int offs_j = atom_j * 3;
      float x1 = coords[offs_j], y1 = coords[offs_j+1], z1 = coords[offs_j+2];
      uint8_t thisBondType = BondType[bondtype_offs_i+atom_j];
      // Skip this atom pair if they are connected by a bond or angle
      if((thisBondType & (BOND|ANGLE)) != 0) continue;

      // Distance reciprocals
      float dx = x1-x0, dy = y1-y0, dz = z1-z0;
      double distRecip = 1 / sqrtf(dx*dx+dy*dy+dz*dz);
      double distRecip3 = distRecip * distRecip * distRecip;
      double distRecip6 = distRecip3 * distRecip3;
      
      // Locate LJ params for this atom pair
      int index = NBIndices[nbparm_offs_i+AtomTypeIndices[atom_j]-1] - 1;
      double thisEnergy = LJ12[index]*distRecip6*distRecip6 - LJ6[index]*distRecip6;

      // Are these atoms 1-4 to each other? If so, divide the energy
			// by 1.2, as ff99 et al dictate.
      if((thisBondType & DIHEDRAL) != 0)
        thisEnergy *= VDW_14_SCALING_RECIP;
      
      decomp[i_res*Nresidues+ResidueMap[atom_j]] += thisEnergy;
      decomp[i_res+ResidueMap[atom_j]*Nresidues] += thisEnergy;
      energy += thisEnergy;
      
    }
  }
  return energy;
  
}

// Calculates pairwise electrostatic energies; no cutoff/PME/etc
double Electro(float *coords, double *decomp) {
  double energy = 0.0;
  int atom_i;
	
  for(atom_i = 0; atom_i < Natoms; atom_i++) {
    int offs_i = atom_i * 3;
    float x0 = coords[offs_i], y0 = coords[offs_i+1], z0 = coords[offs_i+2];
    float qi = Charges[atom_i];
    int i_res = ResidueMap[atom_i]; // Residue of atom i
    // Iterate over all atoms
    int atom_j;
    
    for(atom_j = 0; atom_j < atom_i; atom_j++) {
      int offs_j = atom_j * 3;
      float x1 = coords[offs_j], y1 = coords[offs_j+1], z1 = coords[offs_j+2];
      uint8_t thisBondType = BondType[atom_i*Natoms+atom_j];
      // Skip this atom pair if they are connected by a bond or angle
      if((thisBondType & (BOND|ANGLE)) != 0) continue;

      float dx = x1-x0, dy = y1-y0, dz = z1-z0;
      double thisEnergy = (qi * Charges[atom_j]) / sqrtf(dx*dx+dy*dy+dz*dz);

      // Are these atoms 1-4 to each other? If so, divide the energy
			// by 1.2, as ff99 et al dictate.
      if((thisBondType & DIHEDRAL) != 0)
        thisEnergy *= EEL_14_SCALING_RECIP;
      
      decomp[i_res*Nresidues+ResidueMap[atom_j]] += thisEnergy;
      decomp[i_res+ResidueMap[atom_j]*Nresidues] += thisEnergy;
      energy += thisEnergy;
    }
  }
  return energy;
}

int checksum(char *data, int len) {
  int i;
  int foo = 0;
  for(i = 0; i < len; i++)
    foo += data[i];
  return foo;
}

int main (int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &NumNodes);
  
  // If rank 0, we'll load the data, distribute it to the other nodes, and write stuff to disk
  if(Rank == 0) {
    printf("coulomb (C/MPI version) - T. Joseph <thomas.joseph@mssm.edu>\n\n");
    printf("Running on %d nodes.\n", NumNodes);
    
    if(argc < 4) {
      printf("Usage: <solute.top.tom> <md.trj.gz> <md.ene.bin>\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
    
    FILE *fp = fopen(argv[1], "r");
    if(!fp) bomb("Couldn't open preprocessed prmtop.");
    
    // Load the preprocessed binary version of the prmtop
    fread(&Natoms, sizeof(int), 1, fp);
    fread(&Nresidues, sizeof(int), 1, fp);
    fread(&Ntypes, sizeof(int), 1, fp);
    NBIndices = loadIntArray(fp, &NumNBIndices);
    AtomTypeIndices = loadIntArray(fp, &NumAtomTypeIndices);
    LJ12 = loadFloatArray(fp, &NumLJ12);
    LJ6 = loadFloatArray(fp, &NumLJ6);
    Charges = loadFloatArray(fp, &NumCharges);
    BondType = loadByteArray(fp, &NumBondType);
    ResidueMap = loadIntArray(fp, &NumResidueMap);
    
    if((int)sqrt(NumBondType) != Natoms) {
      printf("sqrt(NumBondType) = %d should be equal to Natoms = %d\n", (int)sqrt(NumBondType), Natoms);
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
    
    // Calculate how many lines per frame in the trajectory file, for ease of reading
    int hasBox = 1; // XXX: We always assume there's a box
    int linesPerFrame = Natoms * 3 / 10;
    if(Natoms*3%10 != 0)
      linesPerFrame++;
    
    // Open the trajectory file
    gzFile trj = gzopen(argv[2], "r");
    char line[256];
    // Eat header line
    gzgets(trj, line, 256);
    
    // Open the md.ene.bin file - for writing raw decomposition matrices
    FILE *eneOut = fopen(argv[3], "w");
  
    // Tell all the worker nodes about the molecule we're examining
    syncMolecule();
    
    // Now that the molecule information is synchronized, we can ask the other nodes to
    // calculate things. Each request will be a float array with (Natoms * 3 + 1) elements.
    // The first float is a command: 0 = Valid request, 1 = Calculation finished so exit.
    // The rest of the floats are the coordinates.
    
    // While there are free nodes, assign calculations to them using MPI_Send and receive
    // results using MPI_Irecv. Use MPI_Waitany to wait until one of the nodes finishes.
    // When it does, write its output to disk.
    int occupied[NumNodes-1];
    memset(occupied, 0, sizeof(int) * (NumNodes-1));
    float *requestBuf[NumNodes-1];
    double *resultBuf[NumNodes-1];
    int node;
    // Allocate request/coordinate buffers
    for(node = 0; node < (NumNodes-1); node++) {
      requestBuf[node] = (float*) malloc(sizeof(float) * (Natoms*3+1));
      requestBuf[node][0] = 0.0f; // Keep on truckin'
      resultBuf[node] = (double*) malloc(sizeof(double) * (Nresidues*Nresidues));
    }
    float *floatBuf = (float*) malloc(sizeof(float) * (Nresidues*Nresidues));
    
    MPI_Request requests[NumNodes-1];
    memset(requests, 0, sizeof(MPI_Request) * (NumNodes-1));
    
    int numFramesDone = 0, finished = 0, occupiedCount = 0;
    
    for(;;) {
      // Assign work to any unassigned nodes
      for(node = 0; node < (NumNodes-1); node++) {
        if(!occupied[node]) {
          // Load a frame from the trajectory file. 10 reals per line, 8 chars per real
          int idx = 1, i, j;
          for(i = 0; i < linesPerFrame; i++) {
            char *ret = gzgets(trj, line, 256);
            if(ret == NULL) {
              // No frames left, so stop
              finished = 1;
              break;
            }
            int len = strlen(line);
            for(j = 0; j < len-1; j+=8) {
              requestBuf[node][idx] = (float) strtod(line+j, NULL);
              idx++;
            }
          }
          if(finished) break; // Stop trying to assign more work
          if(hasBox) gzgets(trj, line, 256); // Eat box info
          occupied[node] = 1; // Mark that node as occupied
          occupiedCount++;
          // requestBuf[node][0] = 0.0f;
          if(MPI_Send(requestBuf[node], Natoms*3+1, MPI_FLOAT, node+1, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
            bomb("MPI_Send assigning a request");
          if(MPI_Irecv(resultBuf[node], Nresidues*Nresidues, MPI_DOUBLE, node+1, 0, MPI_COMM_WORLD, &requests[node]) != MPI_SUCCESS)
            bomb("MPI_Irecv receving results of a request");
        } // if !occupied[node]
        
      } // for(nodes)
    
      int index, waitFor = 1;
      // Collect results from all occupied nodes if we're finished
      if(finished) {
        waitFor = occupiedCount;
        printf("OK, I think we're all finished. Just gotta get those last few frames: %d\n", waitFor);
      }
      MPI_Status status;
      int i;
      for(i = 0; i < waitFor; i++) {
        if(MPI_Waitany(NumNodes-1, requests, &index, &status) != MPI_SUCCESS)
          bomb("MPI_Waitany");
        // node rank 1 is index 0 here
        if(index == MPI_UNDEFINED)
          bomb("MPI_Waitany failed with MPI_UNDEFINED");
          
        numFramesDone++;
        if((numFramesDone % 100) == 0)
          printf("Frames done: %d\n", numFramesDone);
          
        occupied[index] = 0; // Mark that node as unoccupied
        occupiedCount--;
        
        // Convert decomp matrix to floats and dump it to the file
        int k;
        for(k = 0; k < (Nresidues*Nresidues); k++)
          floatBuf[k] = (float) resultBuf[index][k];
        fwrite(floatBuf, sizeof(float), Nresidues*Nresidues, eneOut);
      } // for waitFor
      
      if(finished) {
        // Tell node to quit if we're finished
        for(node = 1; node < NumNodes; node++) {
          requestBuf[node-1][0] = 1.0f;
          if(MPI_Send(requestBuf[node-1], Natoms*3+1, MPI_FLOAT, node, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
            bomb("MPI_Send Telling worker node to quit");
        }
        break; // Quit in general
      }
    }
    fclose(eneOut);
  } else {
    // Worker node stuff
    // Receive prmtop information
    syncMolecule();
    // printf("[%d] Worker node synced and ready to rock!\n", Rank);
    // printf("[%d] Ntypes = %d, Natoms = %d, Nresidues = %d, NumNBIndices = %d, NumCharges = %d, NumBondType = %d\n", Rank, Ntypes, Natoms, Nresidues, NumNBIndices, NumCharges, NumBondType);
    
    // Allocate request buffer
    float *requestBuf = (float*) malloc(sizeof(float) * (Natoms*3+1));
    double *resultBuf = (double*) malloc(sizeof(double) * (Nresidues*Nresidues));
    
    for(;;) {
      // Receive request.
      if(MPI_Recv(requestBuf, Natoms*3+1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
        bomb("Receiving request");
      // printf("[%d] Received request.\n", Rank);
      // If we are done, quit. Else, process and send results back; repeat.
      if(requestBuf[0] != 0.0f) {
        printf("[%d] I was told to quit, so see ya!\n", Rank);
        break;
      }
      float *coords = requestBuf+1;

      memset(resultBuf, 0, sizeof(double) * (Nresidues*Nresidues));
      double eel = Electro(coords, resultBuf);
      double vdw = LennardJones(coords, resultBuf);
      
      // DEBUG
      int r, c;
      double ene = 0.0;
      for(r = 0; r < Nresidues; r++) {
        for(c = 0; c <= r; c++) {
          ene += (double) resultBuf[r*Nresidues+c];
        }
      }
      
      printf("[%d] Electro: %f vdW: %f Total: %f FakeTotal: %f\n", Rank, eel, vdw, eel+vdw, ene);
      if(MPI_Send(resultBuf, Nresidues*Nresidues, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
        bomb("Sending back results");
    }
  }
  
  MPI_Finalize();  
  return 0;
}