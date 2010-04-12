#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <math.h>
#include <mpi.h>

#define CUTOFF 10.0

int Rank, NumNodes;
int Nresidues;

// Synchronizes data by broadcasting it from the rank 0 node
void syncIntArray(int **buf, int *size) {
  MPI_Bcast(size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  if(Rank > 0) // If we're not the master node, we need to allocate this buffer
    *buf = malloc(sizeof(int) * *size);
  MPI_Bcast(*buf, *size, MPI_INTEGER, 0, MPI_COMM_WORLD);
}

int numFramesInFile(char *filename, int numResidues) {
  int frameSize = numResidues*numResidues*sizeof(float);
  struct stat info;
  stat(filename, &info);
  return (int) (info.st_size / frameSize);
}

int main (int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &NumNodes);
  
  // If rank 0, we'll load the data, distribute it to the other nodes, and write stuff to disk
  if(Rank == 0) {
    printf("correl (C/MPI version) - T. Joseph <thomas.joseph@mssm.edu>\n\n");
    printf("Running on %d nodes.\n", NumNodes);
    
    // Usage: <num-residues> <md.ene.bin.0> [md.ene.bin.1] [...]
    if(argc < 3) {
      printf("Usage: <num-residues> <md.ene.bin.0> [md.ene.bin.1] [...]\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
    
    Nresidues = atoi(argv[1]);
    printf("You said there are %d residues. I hope that's correct.\n", Nresidues);
    
    // Calculate average energies
    // Create and zero out the average energies matrix
    double *averageEnergies = malloc(Nresidues*Nresidues*sizeof(double));
    memset(averageEnergies, 0, Nresidues*Nresidues*sizeof(double));
    
    int i;
    int totalNumFrames = 0;
    float *buf = malloc(Nresidues*Nresidues*sizeof(float));
    for(i = 2; i < argc; i++) {
      int numFrames = numFramesInFile(argv[i], Nresidues);
      printf("Calculating average matrix: %s with %d frames...\n", argv[i], numFrames);
      FILE *fp = fopen(argv[i], "r");
      int frame;
      for(frame = 0; frame < numFrames; frame++) {
        // Load a matrix and add it to the running total
        fread(buf, sizeof(float), Nresidues*Nresidues, fp);
        int j;
        for(j = 0; j < Nresidues*Nresidues; j++)
          averageEnergies[j] += (double) buf[j];
      }
      fclose(fp);
      totalNumFrames += numFrames;
    }
    
    // Divide matrix by number of frames and then count the residue pairs above the cutoff.
    // Ideally we'd use some sort of resizing vector to store them, but instead we'll use the gimpy
    // technique of counting how many there are on a first pass, then allocating an array of the
    // correct size, then actually storing the indices on a second pass. The matrix is small enough
    // that doubling the runtime won't be noticed on today's monster computers.
    int j, numPairs = 0;
    for(j = 0; j < Nresidues*Nresidues; j++) {
      averageEnergies[j] /= (double) totalNumFrames;
    }
    
    int res_i, res_j;
    for(res_i = 0; res_i < Nresidues; res_i++) {
      for(res_j = 0; res_j < res_i; res_j++) {
        if(fabs(averageEnergies[Nresidues*res_i+res_j]) >= CUTOFF)
          numPairs++;
      }
    }
    
    printf("Found %d pairs above cutoff of %.1f kcal/mol.\n", numPairs, CUTOFF);
    
    int *pairsList = malloc(2*numPairs*sizeof(int));
    int ptr = 0;
    for(res_i = 0; res_i < Nresidues; res_i++) {
      for(res_j = 0; res_j < res_i; res_j++) {
        if(averageEnergies[Nresidues*res_i+res_j] >= CUTOFF) {
          pairsList[ptr] = res_i;
          pairsList[ptr+1] = res_j;
          ptr+=2;
        }
      }
    }
    
    // OK! Now we have the pairs list. Now to calculate the correlation matrix (numPairs*numPairs).
    // Broadcast the pairs list.
    int pairsListSize = 2*numPairs;
    syncIntArray(&pairsList, &pairsListSize);
    
    double *num = malloc(Nresidues*Nresidues*sizeof(double));
    double *denom = malloc(Nresidues*Nresidues*sizeof(double));
    memset(num, 0, Nresidues*Nresidues*sizeof(double));
    memset(denom, 0, Nresidues*Nresidues*sizeof(double));
    
    // Load a batch of frames from the current file. We can't just inhale an entire file 
    // because we probably don't have enough RAM to do so. If we're out of frames in the
    // current file, open a new one. If we're out of files, move on to dividing num by denom.
    
    // Make the raw data for correlation calculation, as follows:
    // pair0 ene0 ene1 ene2 ...
    // pair1 ene0 ene1 ene2 ...
    // pair2 ene0 ene1 ene2 ...
    //  ...
    // That is, numPairs rows and numFrames columns in a matrix
    int pair, frame;
    for(pair = 0; pair < numPairs; pair++) {
      int pair_i = pairsList[pair*2];
      int pair_j = pairsList[pair*2+1];
      // Iterate through frames
    }
    
    // Broadcast this pairs-energies matrix.
    
    // Each worker node gets to calculate a roughly equal share of rows of the correlation matrix.
    // The last node gets the modulus of rows left over from the integer division.

    
  } else {
    // Worker node

    // Broadcast the pairs list from master.
    int *pairsList, pairsListSize;
    syncIntArray(&pairsList, &pairsListSize);
    
    // Repeat:
    // Broadcast the pairs-energies matrix from master.
    
    // Receive the range of rows to calculate, and do so. Or else a signal to move on to the next phase.
    
    // Send back the rows.
    // End repeat.
    
    // Broadcast the finished correlation matrix from master.
    
    // Receive range of residues to do the projection with (e.g., range of rows).
    // Project the correlation matrix back onto the residues, resulting in residue interaction matrix.
    // Send back the rows of the residue interaction matrix we just calculated.

    // Finished!
  }
  MPI_Finalize();
  return 0;
}