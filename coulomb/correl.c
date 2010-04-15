#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <math.h>
#include <mpi.h>

#define CUTOFF 10.0
#define FRAME_BATCH_SIZE 500

int Rank, NumNodes;
int Nresidues;

void bomb(char *msg) {
  printf("[%d] ERROR, aborting! %s\n", Rank, msg);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}

// Synchronizes data by broadcasting it from the rank 0 node
void syncIntArray(int **buf, int *size) {
  MPI_Bcast(size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  if(Rank > 0) // If we're not the master node, we need to allocate this buffer
    *buf = malloc(sizeof(int) * *size);
  MPI_Bcast(*buf, *size, MPI_INTEGER, 0, MPI_COMM_WORLD);
}

size_t numFramesInFile(char *filename, int numResidues) {
  int frameSize = numResidues*numResidues*sizeof(float);
  struct stat info;
  stat(filename, &info);
  return (size_t) (info.st_size / frameSize);
}

int min(int a, int b) {
  if(a < b)
    return a;
  return b;
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
    MPI_Bcast(&Nresidues, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    
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

    // Workers need this to calculate correlation matrix
    MPI_Bcast(averageEnergies, Nresidues*Nresidues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
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
    MPI_Bcast(rowsToCalculate, NumNodes, MPI_INTEGER, 0, MPI_COMM_WORLD);
    
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
        printf("Calculating correlation matrix: On file %s, frames %d to %d.\n", argv[fileIdx], i, i + thisBlockSize);
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
            pairsEnergies[thisBlockSize * pair + frame] = buf[numPairs*pair_i+pair_j];
          }
        } // iterating over frames within block
        
        // At this point, we have a complete pairsEnergies matrix from which we can figure some of
        // the correlation matrix. That's a job for the worker nodes.
        // There are a zillion MPI calls below but they are executed pretty infrequently, so I am
        // not bothering with coalescing them.
        keepGoing = 1; // Keep on truckin' - we've got frames to process
        MPI_Bcast(&keepGoing, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
        // Broadcast this pairs-energies matrix. First, how many frames: matrix width.
        // Matrix height is numPairs. Then, the data itself.
        MPI_Bcast(&thisBlockSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
        MPI_Bcast(pairsEnergies, thisBlockSize*numPairs, MPI_FLOAT, 0, MPI_COMM_WORLD);
        
        
      } // iterating over frame blocks

      // Done with this file
      fclose(fp);
    }

    // Done with all files
    keepGoing = 0;
    MPI_Bcast(&keepGoing, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    // Free up some memory
    free(pairsEnergies);
    
    // Receive numerator and denominator matrix parts that the worker nodes calculated.
    // We do this serially because
    //   a) I couldn't get MPI_Irecv/MPI_Waitall to work properly
    //   b) Latency is not important at this stage anyhow
    double *numLocal[NumNodes];
    double *denomLocal[NumNodes];
    MPI_Request numRequests[NumNodes], denomRequests[NumNodes];
    memset(numRequests, 0, sizeof(MPI_Request) * NumNodes);
    memset(denomRequests, 0, sizeof(MPI_Request) * NumNodes);
    for(node = 1; node < NumNodes; node++) {
      numLocal[node] = malloc(sizeof(double) * localSizes[node]);
      denomLocal[node] = malloc(sizeof(double) * localSizes[node]);
      MPI_Recv(numLocal[node], localSizes[node], MPI_DOUBLE, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(denomLocal[node], localSizes[node], MPI_DOUBLE, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Assemble num and denom from the pieces we just received
    double *num = malloc(numPairs*numPairs*sizeof(double));
    double *denom = malloc(numPairs*numPairs*sizeof(double));
    // ...
    
    // Projection! Sometime in the future probably! For now the Go implementation is
    // probably fast enough.
    
  } else {
    // We are a worker node. Our mission is to serve.
    MPI_Bcast(&Nresidues, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    double *averageEnergies = malloc(Nresidues*Nresidues*sizeof(double));
    MPI_Bcast(averageEnergies, Nresidues*Nresidues, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Broadcast the pairs list from master.
    int *pairsList, pairsListSize;
    syncIntArray(&pairsList, &pairsListSize);
    int numPairs = pairsListSize/2;
    // printf("[%d] Hi there. I was told there are %d residues and %d pairs.\n", Rank, Nresidues, numPairs);
    float *pairsEnergies = malloc(numPairs*FRAME_BATCH_SIZE*sizeof(float));
    int keepGoing;

    // Receive the range of rows to calculate (ij_a, ij_b).
    int *rowsToCalculate = malloc(sizeof(int) * NumNodes);
    MPI_Bcast(rowsToCalculate, NumNodes, MPI_INTEGER, 0, MPI_COMM_WORLD);
    int ij_a = rowsToCalculate[Rank-1], ij_b = rowsToCalculate[Rank];
    // Easy way - this wastes space - see below
    int localSize = (ij_b-ij_a) * ij_b;
    printf("[%d] ij_a = %d, ij_b = %d, localSize = %.0fMB.\n", Rank, ij_a, ij_b, (double)localSize*sizeof(double)/1048576.0);
    double *numLocal = malloc(localSize * sizeof(double));
    double *denomLocal = malloc(localSize * sizeof(double));
    memset(numLocal, 0, localSize * sizeof(double));
    memset(denomLocal, 0, localSize * sizeof(double));
    
    // Repeat:
    for(;;) {
      // Should we keep going?
      MPI_Bcast(&keepGoing, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      if(!keepGoing) break;

      // Broadcast the pairs-energies matrix from master.
      int thisBlockSize;
      MPI_Bcast(&thisBlockSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      MPI_Bcast(pairsEnergies, thisBlockSize*numPairs, MPI_FLOAT, 0, MPI_COMM_WORLD);
      // printf("[%d] Received pairsEnergies matrix with %d floats in it.\n", Rank, thisBlockSize*numPairs);
    
      // Actually do the calculation.
      int ij, kl;
      for(ij = ij_a; ij < ij_b; ij++) {
        for(kl = 0; kl < ij; kl++) {
          int i = pairsList[ij*2], j = pairsList[ij*2+1];
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
            //   M[numPairs*ij+kl] -> M'[(ij_a+m)*m+kl] where m = ij-ij_a
            // This indexing method wastes space, but whatever, it's simple. We take all
            // rows to have width ij_b and just don't use the upper half on the right side
            // of our local slab. We could do fancier indexing to save that space but the
            // user won't notice all the extra effort we went to.
            int m = ij-ij_a; // "local" row index
            int offs = ij_b*m+kl;
            numLocal[offs] += a*b;
            denomLocal[offs] += sqrt(a*a*b*b);
          }
        }
      }
    
      // End repeat.
    }

    // Don't need this anymore
    free(averageEnergies);
    // Send back our local piece of num and denom
    MPI_Send(numLocal, localSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free(numLocal); // Free up space in case the master node is running on this machine
    MPI_Send(denomLocal, localSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free(denomLocal);
      
    // For the future, potentially:
    //   Broadcast the finished correlation matrix from master. (That could be a problem because
    //     it could be nearly 2GB in size.)
    //   Receive range of residues to do the projection with (e.g., range of rows).
    //   Project the correlation matrix back onto the residues, resulting in residue interaction matrix.
    //   Send back the rows of the residue interaction matrix we just calculated.

    // Finished!
  }
  MPI_Finalize();
  return 0;
}