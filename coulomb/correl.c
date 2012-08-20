#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#include <math.h>
#include <mpi.h>

#define ENERGY_CUTOFF 10.0
#define CORREL_CUTOFF 0.4
#define FRAME_BATCH_SIZE 500
#define MAX_HOSTNAME_LENGTH 1024

char Hostname[MAX_HOSTNAME_LENGTH];
int Rank, NumNodes;
int Nresidues;

void bomb(char *msg) {
  printf("[%d:%s] ERROR, aborting! %s\n", Rank, Hostname, msg);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}

// Synchronizes data by broadcasting it from the rank 0 node
void broadcastIntArray(int **buf, int *size) {
  MPI_Bcast(size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(Rank > 0) // If we're not the master node, we need to allocate this buffer
    *buf = malloc(sizeof(int) * *size);
  MPI_Bcast(*buf, *size, MPI_INT, 0, MPI_COMM_WORLD);
}

// Returns the number of frames in a file with numResidues*numResidues*sizeof(float) bytes per frame
size_t numFramesInFile(char *filename, int numResidues) {
  int frameSize = numResidues*numResidues*sizeof(float);
  struct stat info;
  stat(filename, &info);
  return (size_t) (info.st_size / frameSize);
}

// Returns minimum of integers a and b
int min(int a, int b) {
  if(a < b)
    return a;
  return b;
}

// Assumes matrix is stored row by row
void dumpIntMatrixToFile(char *filename, int *data, int numRows, int numCols) {
  FILE *fp = fopen(filename, "w");
  int row, col, ptr = 0;
  for(row = 0; row < numRows; row++) {
    for(col = 0; col < numCols; col++) {
      if(col != 0) fprintf(fp, " ");
      fprintf(fp, "%d", data[ptr++]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

// Dumps a matrix of doubles to a text file
void dumpDoubleMatrixToFile(char *filename, double *data, int numRows, int numCols) {
  FILE *fp = fopen(filename, "w");
  int row, col, ptr = 0;
  for(row = 0; row < numRows; row++) {
    for(col = 0; col < numCols; col++) {
      if(col != 0) fprintf(fp, " ");
      fprintf(fp, "%f", data[ptr++]);
    }
    fprintf(fp, "\n");
  }
  fflush(fp);
  fclose(fp);
}

// Dumps a matrix of floats to a text file
void dumpFloatMatrixToFile(char *filename, float *data, int numRows, int numCols) {
  FILE *fp = fopen(filename, "w");
  int row, col, ptr = 0;
  for(row = 0; row < numRows; row++) {
    for(col = 0; col < numCols; col++) {
      if(col != 0) fprintf(fp, " ");
      fprintf(fp, "%f", data[ptr++]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

int main (int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  gethostname(Hostname, MAX_HOSTNAME_LENGTH);
  MPI_Comm_size(MPI_COMM_WORLD, &NumNodes);
  
  // If rank 0, we'll load the data, distribute it to the other nodes, and write stuff to disk
  if(Rank == 0) {
    printf("correl (C/MPI version) on %d nodes - T. Joseph <thomas.joseph@mssm.edu>\n\n", NumNodes);
    printf("Running on %d nodes.\n", NumNodes);
    
    // Usage: <num-residues> <md.ene.bin.0> [md.ene.bin.1] [...]
    if(argc < 3) {
      printf("Usage: <num-residues> <md.ene.bin.0> [md.ene.bin.1] [...]\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
    
    Nresidues = atoi(argv[1]);
    printf("You said there are %d residues. I hope that's correct.\nUsing energy cutoff of %.1f kcal/mol and correlation cutoff of %.1f.\n", Nresidues, ENERGY_CUTOFF, CORREL_CUTOFF);
    if(MPI_Bcast(&Nresidues, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      bomb("Broadcasting number of residues");
      
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
        int r, c;
        double ene = 0.0;
        for(r = 0; r < Nresidues; r++) {
          for(c = 0; c <= r; c++) {
            ene += (double) buf[r*Nresidues+c];
          }
        }
        // printf("Frame %d energy: %f\n", frame, ene); // DEBUG
        
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
    int j;
    for(j = 0; j < Nresidues*Nresidues; j++) {
      averageEnergies[j] /= (double) totalNumFrames;
    }
    
    dumpDoubleMatrixToFile("average.txt", averageEnergies, Nresidues, Nresidues);
    printf("Dumped average matrix to average.txt.\n");

    // Workers need this to calculate correlation matrix
    if(MPI_Bcast(averageEnergies, Nresidues*Nresidues, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      bomb("Broadcasting averageEnergies matrix");
    
    int res_i, res_j, numPairs = 0;
    for(res_i = 0; res_i < Nresidues; res_i++) {
      for(res_j = 0; res_j < (res_i-1); res_j++) {
        if(fabs(averageEnergies[Nresidues*res_i+res_j]) >= ENERGY_CUTOFF)
          numPairs++;
      }
    }
    
    int *pairsList = malloc(2*numPairs*sizeof(int)), ptr = 0;
    for(res_i = 0; res_i < Nresidues; res_i++) {
      for(res_j = 0; res_j < (res_i-1); res_j++) {
        if(fabs(averageEnergies[Nresidues*res_i+res_j]) >= ENERGY_CUTOFF) {
          pairsList[ptr] = res_i;
          pairsList[ptr+1] = res_j;
          ptr+=2;
        }
      }
    }
    
    printf("Found %d pairs above cutoff of %.1f kcal/mol.\n", numPairs, ENERGY_CUTOFF);
    
    // OK! Now we have the pairs list. Now to calculate the correlation matrix (numPairs*numPairs).
    
    // Broadcast the pairs list.
    int pairsListSize = 2*numPairs;
    broadcastIntArray(&pairsList, &pairsListSize);
    
    float *pairsEnergies = malloc(numPairs*FRAME_BATCH_SIZE*sizeof(float));
    // Once we know numPairs we can determine rowsToCalculate, so we only have to send it once.
    // Each worker node gets to calculate a roughly equal share of the correlation matrix.
    int *rowsToCalculate = malloc(NumNodes * sizeof(int));
    int *localSizes = malloc(NumNodes * sizeof(int));
    rowsToCalculate[0] = 0; // Start calculating at row 0
    int node;
    // Can't just naively divide - have to take into account the shape of the share we're assigning.    
    for(node = 1; node < NumNodes; node++) {
      rowsToCalculate[node] = (int) sqrt(((double)node * (double)numPairs * (double)numPairs) / (double)(NumNodes-1));
      // int localSize = (ij_b-ij_a) * ij_b;
      localSizes[node] = (rowsToCalculate[node] - rowsToCalculate[node-1]) * rowsToCalculate[node];
    }
    // Make sure rounding error doesn't leave a row uncalculated
    rowsToCalculate[NumNodes] = numPairs;
    if(MPI_Bcast(rowsToCalculate, NumNodes, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      bomb("Broadcasting rowsToCalculate");
    
    int fileIdx, keepGoing;
    for(fileIdx = 2; fileIdx < argc; fileIdx++) {
      size_t numFrames = numFramesInFile(argv[fileIdx], Nresidues);
      // numFrames = 50; // DEBUG
      FILE *fp = fopen(argv[fileIdx], "r");
      // Load a batch of frames from the current file. We can't just inhale an entire file 
      // because we probably don't have enough RAM to do so.
      int i;
      for(i = 0; i < numFrames; i += FRAME_BATCH_SIZE) {
        int thisBlockSize = min(numFrames - i, FRAME_BATCH_SIZE);
        int frame;
        printf("Calculating pairs correlation matrix: On file %s, frames %d to %d.\n", argv[fileIdx], i, i + thisBlockSize);
        for(frame = 0; frame < thisBlockSize; frame++) {
          // Load a frame
          fread(buf, sizeof(float), Nresidues*Nresidues, fp);
          
          // Make the raw data for correlation calculation, as follows:
          // pair0: ene0 ene1 ene2 ...
          // pair1: ene0 ene1 ene2 ...
          // pair2: ene0 ene1 ene2 ...
          //  ...
          // That is, numPairs rows and numFrames columns in a matrix
          int pair;
          for(pair = 0; pair < numPairs; pair++) {
            int pair_i = pairsList[pair*2];
            int pair_j = pairsList[pair*2+1];
            pairsEnergies[thisBlockSize * pair + frame] = buf[Nresidues*pair_i+pair_j];
          }
        } // iterating over frames within block
        
        // At this point, we have a complete pairsEnergies matrix from which we can figure some of
        // the correlation matrix. That's a job for the worker nodes.
        // There are a zillion MPI calls below but they are executed pretty infrequently, so I am
        // not bothering with coalescing them.
        keepGoing = 1; // Keep on truckin' - we've got frames to process
        if(MPI_Bcast(&keepGoing, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
          bomb("Broadcasting keepGoing");
        // Broadcast this pairs-energies matrix. First, how many frames: matrix width.
        // Matrix height is numPairs. Then, the data itself.
        if(MPI_Bcast(&thisBlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
          bomb("Broadcasting thisBlockSize");
        if(MPI_Bcast(pairsEnergies, thisBlockSize*numPairs, MPI_FLOAT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
          bomb("Broadcasting pairsEnergies");
      } // iterating over frame blocks

      // Done with this file
      fclose(fp);
    }

    // Done with all files
    keepGoing = 0;
    if(MPI_Bcast(&keepGoing, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      bomb("Broadcasting stop signal");
    // Free up some memory
    free(pairsEnergies);
    
    // Receive numerator and denominator matrix parts that the worker nodes calculated.
    // We do this serially because
    //   a) I couldn't get MPI_Irecv/MPI_Waitall to work properly
    //   b) Latency is not important at this stage anyhow
    // Also: Assemble num and denom from the pieces we just received
    double *num = malloc(numPairs*numPairs*sizeof(double));
    double *denom = malloc(numPairs*numPairs*sizeof(double));
    for(node = 1; node < NumNodes; node++) {
      double *numLocal = malloc(sizeof(double) * localSizes[node]);
      double *denomLocal = malloc(sizeof(double) * localSizes[node]);
      MPI_Recv(numLocal, localSizes[node], MPI_DOUBLE, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(denomLocal, localSizes[node], MPI_DOUBLE, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // M[numPairs*ij+kl] <- M'[ij_b*m+kl] where m = ij-ij_a
      size_t ij, kl, ij_a = rowsToCalculate[node-1], ij_b = rowsToCalculate[node];
      for(ij = ij_a; ij < ij_b; ij++) {
        size_t m = ij-ij_a; // "local" row index
        for(kl = 0; kl < ij; kl++) {
          size_t offs = ij_b*m+kl;
          num[numPairs*ij+kl] = num[numPairs*kl+ij] = numLocal[offs];
          denom[numPairs*ij+kl] = denom[numPairs*kl+ij] = denomLocal[offs];
        }
      }
      free(numLocal);
      free(denomLocal);
    }
    
    // Now that we've reassembled num and denom, divide them to give final correlation matrix.
    // This can be floats because we're just going to dump it to a file, and we want to save memory.
    float *correl = (float*) malloc(numPairs*numPairs*sizeof(float));
    for(i = 0; i < (numPairs*numPairs); i++)
      correl[i] = (float)(num[i]/denom[i]);
      
    // Apply correlation cutoff as in Kong and Karplus.
    // All correlation values that do not meet the cutoff are taken to be zero.
    // Rows/columns that are entirely zero are discarded from the matrix, and the corresponding
    // pairs from pairsList.
    int *savedPairIDs = malloc(numPairs*sizeof(int)), prunedCount = 0;
    for(int ij = 0; ij < numPairs; ij++) {
      for(int kl = 0; kl < numPairs; kl++) {
        // If there is a correl value above the cutoff, this row/column is good so
        // save this index and move on to the next row immediately
        if(ij != kl && fabs(correl[numPairs*ij+kl]) > CORREL_CUTOFF) {
          savedPairIDs[prunedCount++] = ij;
          break;
        }
      }
    }
    
    // Extract the relevant pairs
    int *prunedPairsList = (int*) malloc(prunedCount*2*sizeof(int));
    for(int i = 0; i < prunedCount; i++) {
      prunedPairsList[i*2] = pairsList[savedPairIDs[i]*2];
      prunedPairsList[i*2+1] = pairsList[savedPairIDs[i]*2+1];
    }
    
    // Extract the relevant correlation values and make a new matrix
    float *prunedCorrel = malloc(prunedCount*prunedCount*sizeof(float));
    for(int ij = 0; ij < prunedCount; ij++) {
      for(int kl = 0; kl < prunedCount; kl++) {
        // savedPairIDs contains the indices into the original correl matrix
        int orig_ij = savedPairIDs[ij], orig_kl = savedPairIDs[kl];
        prunedCorrel[prunedCount*ij+kl] = correl[numPairs*orig_ij+orig_kl];
      }
    }
      
    dumpIntMatrixToFile("orig-pairs.txt", pairsList, numPairs, 2);
    printf("Dumped unpruned list of residue pairs to orig-pairs.txt.\n");
    dumpFloatMatrixToFile("orig-correl.txt", correl, numPairs, numPairs);
    printf("Dumped unpruned pairs correlation matrix to orig-correl.txt.\n");

    dumpIntMatrixToFile("pairs.txt", prunedPairsList, prunedCount, 2);
    printf("Dumped pruned list of residue pairs to pairs.txt (those that met the correlation cutoff of %.1f).\n", CORREL_CUTOFF);
    dumpFloatMatrixToFile("correl.txt", prunedCorrel, prunedCount, prunedCount);
    printf("Dumped pruned pairs correlation matrix to correl.txt.\n");
    
  } else {
    // We are a worker node. Our mission is to serve.
    MPI_Bcast(&Nresidues, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double *averageEnergies = malloc(Nresidues*Nresidues*sizeof(double));
    MPI_Bcast(averageEnergies, Nresidues*Nresidues, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Broadcast the pairs list from master.
    int *pairsList, pairsListSize;
    broadcastIntArray(&pairsList, &pairsListSize);
    int numPairs = pairsListSize/2;
    // printf("[%d] Hi there. I was told there are %d residues and %d pairs.\n", Rank, Nresidues, numPairs);
    float *pairsEnergies = malloc(numPairs*FRAME_BATCH_SIZE*sizeof(float));
    int keepGoing;

    // Receive the range of rows to calculate (ij_a, ij_b).
    int *rowsToCalculate = malloc(sizeof(int) * NumNodes);
    MPI_Bcast(rowsToCalculate, NumNodes, MPI_INT, 0, MPI_COMM_WORLD);
    int ij_a = rowsToCalculate[Rank-1], ij_b = rowsToCalculate[Rank];
    // Easy way - this wastes space - see below
    int localSize = (ij_b-ij_a) * ij_b;
    printf("[%d:%s] ij_a = %d, ij_b = %d, localSize = %.0fMB.\n", Rank, Hostname, ij_a, ij_b, (double)localSize*sizeof(double)/1048576.0);
    double *numLocal = malloc(localSize * sizeof(double));
    double *denomLocal = malloc(localSize * sizeof(double));
    memset(numLocal, 0, localSize * sizeof(double));
    memset(denomLocal, 0, localSize * sizeof(double));
    
    // Repeat:
    for(;;) {
      // Should we keep going?
      MPI_Bcast(&keepGoing, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(!keepGoing) break;

      // Broadcast the pairs-energies matrix from master.
      int thisBlockSize;
      MPI_Bcast(&thisBlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(pairsEnergies, thisBlockSize*numPairs, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
      // Actually do the calculation.
      int ij, kl;
      for(ij = ij_a; ij < ij_b; ij++) {
        int m = ij-ij_a; // "local" row index
        int i = pairsList[ij*2], j = pairsList[ij*2+1];
        for(kl = 0; kl < ij; kl++) {
          int k = pairsList[kl*2], l = pairsList[kl*2+1];
          int t;
          for(t = 0; t < thisBlockSize; t++) {
            double a = (double) pairsEnergies[thisBlockSize*ij+t] - averageEnergies[Nresidues*i+j];
            double b = (double) pairsEnergies[thisBlockSize*kl+t] - averageEnergies[Nresidues*k+l];
            // If we had the entire num and denom matrices, we could just do
            //   num[numPairs*ij+kl] += a*b;
            //   denom[numPairs*ij+kl] += sqrt(a*a*b*b);
            // But we don't have num and denom, so we need to convert the location (ij, kl) into an
            // index into our local piece of the num and denom matrices and store these values there.
            //   M[numPairs*ij+kl] -> M'[ij_b*m+kl] where m = ij-ij_a
            // This indexing method wastes space, but whatever, it's simple. We take all
            // rows to have width ij_b and just don't use the upper half on the right side
            // of our local slab. We could do fancier indexing to save that space but the
            // user won't notice all the extra effort we went to.
            int offs = ij_b*m+kl;
            numLocal[offs] += a*b;
            denomLocal[offs] += sqrt(a*a*b*b);
          }
        }
        // Really verbose debugging message
        // printf("[%d:%s] Now on %d, going from %d to %d.\n", Rank, Hostname, ij, ij_a, ij_b);
      }
    
      // End repeat.
    }

    // Don't need these anymore
    free(averageEnergies);
    free(pairsEnergies);
    // Send back our local piece of num and denom
    MPI_Send(numLocal, localSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free(numLocal); // Free up space in case the master node is running on this machine
    MPI_Send(denomLocal, localSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free(denomLocal);
      
    // Finished!
  }
  MPI_Finalize();
  return 0;
}
