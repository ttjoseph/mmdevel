#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include <mpi.h>

#define LINE_BUF_SIZE (2*1024*1024)

int Rank, NumNodes; // For MPI

int NumPairs, NumResidues, *Pairs;
float *Correl, *Rescorrel;


void bomb(char *msg) {
  printf("[%d] ERROR, aborting! %s\n", Rank, msg);
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

void broadcastFloatArray(float **buf, int *size) {
  MPI_Bcast(size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(Rank > 0) // If we're not the master node, we need to allocate this buffer
    *buf = malloc(sizeof(float) * *size);
  MPI_Bcast(*buf, *size, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

// Does 1/NUM_THREADS part of the work of projecting the correlation matrix back onto a
// residue-residue matrix.
void* projectPiece() {
  int res_i = 0, res_j = 0;
  for(res_i = Rank; res_i < NumResidues; res_i += NumNodes) {
    for(res_j = 0; res_j < res_i; res_j++) {
      // Like projectCorrelSingleResiduePair
      double totalCorr = 0.0;
      int ij, kl;
      for(ij = 0; ij < NumPairs; ij++) {
        int pair_i = Pairs[ij*2];
        int pair_j = Pairs[ij*2+1];
        for(kl = 0; kl < ij; kl++) {
          int pair_k = Pairs[kl*2];
          int pair_l = Pairs[kl*2+1];
          if(((pair_i == res_i || pair_j == res_i) && (pair_k == res_j || pair_l == res_j))
            || ((pair_i == res_j || pair_j == res_j) && (pair_k == res_i || pair_l == res_i))) {
              totalCorr += (double) fabsf(Correl[NumPairs*ij+kl]);
            }      
        } // ij
      } // kl
      
      Rescorrel[NumResidues*res_i+res_j] = Rescorrel[NumResidues*res_j+res_i] = (float) totalCorr;
    } // res_j
    
    // Estimate how much we've done, assuming that the size of each row grows by one each time
    float percent = (float)res_i;
    percent *= percent;
    percent /= (float) (NumResidues*NumResidues);
    printf("[%d] Finished row %d of %d, about %.0f%% done.\n", Rank, res_i+1, NumResidues, percent*100);
  } // res_i
  
  return NULL;
}

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

int Master(int argc, char *argv[]) {
  printf("projection (C/MPI version) - T. Joseph <thomas.joseph@mssm.edu>\n\n");
  
  if(argc < 2) {
    bomb("Usage: <num-residues>\n\nWe assume you want to load correl.txt and pairs.txt and write to rescorrel.txt.\n");
  }

  // Open and read the first line of correl.txt to decide on its dimensions
  FILE *fp = fopen("correl.txt", "r");
  char *buf = malloc(LINE_BUF_SIZE);
  fgets(buf, LINE_BUF_SIZE, fp);
  
  NumPairs = 1;
  strtok(buf, " "); // Assume there's at least one token
  while(strtok(NULL, " "))
    NumPairs++;
  
  printf("I think there are %d pairs represented in correl.txt, for a structure with %d residues.\n", NumPairs,
    NumResidues);
  rewind(fp);
  
  // Now we know how big the matrix is, so we can inhale it and begin processing
  printf("Allocating %.0fMB for the correl matrix. Loading...\n", ((double) NumPairs)*NumPairs*sizeof(float)/1048576);
  Correl = malloc(sizeof(float)*NumPairs*NumPairs);
  Pairs = malloc(sizeof(int)*NumPairs*2);
  
  // Really not too robust file reading...assumes perfectly formed correl.txt
  int line, tok;
  size_t idx = 0;
  for(line = 0; line < NumPairs; line++) {
    fgets(buf, LINE_BUF_SIZE, fp);
    Correl[idx++] = (float) strtod(strtok(buf, " "), NULL);
    for(tok = 1; tok < NumPairs; tok++) {
      Correl[idx++] = (float) strtod(strtok(NULL, " "), NULL);
    }
  }
  fclose(fp);
  if(idx != (NumPairs*NumPairs)) {
    bomb("ERROR: We didn't load as many numbers from correl.txt as we should have.\n");
  }
  
  // Pretty ghetto to copy and paste code like this, but oh well
  fp = fopen("pairs.txt", "r");
  idx = 0;
  for(line = 0; line < NumPairs; line++) {
    fgets(buf, LINE_BUF_SIZE, fp);
    Pairs[idx++] = atoi(strtok(buf, " "));
    Pairs[idx++] = atoi(strtok(NULL, " "));
    // printf("%d %d\n", Pairs[idx-2], Pairs[idx-1]); // DEBUG
  }
  fclose(fp);
  free(buf);
  
  // Broadcast Correl and Pairs
  int correlMatrixSize = NumPairs*NumPairs;
  int pairsArraySize = NumPairs*2;
  broadcastFloatArray(&Correl, &correlMatrixSize);
  broadcastIntArray(&Pairs, &pairsArraySize);
  
  // As the master, we don't need the Correl matrix and Pairs anymore - just Rescorrel
  free(Correl);
  free(Pairs);
  // We also need a buffer to receive each worker's Rescorrel
  float *RescorrelFromWorker = malloc(sizeof(float)*NumResidues*NumResidues);
  bzero(RescorrelFromWorker, sizeof(float)*NumResidues*NumResidues);
  
  // Receive and combine the parts of Rescorrel from the workers
  for(int workerRank=1; workerRank<NumNodes; workerRank++) {
    MPI_Status status;
    if(MPI_Recv(RescorrelFromWorker, NumResidues*NumResidues, MPI_FLOAT, workerRank, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
      bomb("Receiving result from worker");
    
    for(int i=0; i<(NumResidues*NumResidues); i++)
      Rescorrel[i] += RescorrelFromWorker[i];
    
    printf("Received and combined from worker %d", workerRank);
  }
    
  dumpFloatMatrixToFile("rescorrel.txt", Rescorrel, NumResidues, NumResidues);
  printf("Saved residue correlation matrix to rescorrel.txt. All done, enjoy.\n");
  
  return 0;
}

int Worker(int argc, char *argv[]) {
  // Receive Correl and Pairs
  int correlMatrixSize = NumPairs*NumPairs;
  int pairsArraySize = NumPairs*2;
  broadcastFloatArray(&Correl, &correlMatrixSize);
  broadcastIntArray(&Pairs, &pairsArraySize);
  
  NumPairs = pairsArraySize / 2;
  printf("[%d] Got the goods for %d residues and %d pairs. Let's go!\n", Rank, NumResidues, NumPairs);
  
  projectPiece();
  
  // Send back results to the master process
  if(MPI_Send(Rescorrel, NumResidues*NumResidues, MPI_FLOAT, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    bomb("Sending back results");
  
  return 0;
}

int main (int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &NumNodes);

  // Each process will have its own Rescorrel matrix. The worker matrices will be partial and the
  // master will combine them all to yield the final result
  NumResidues = atoi(argv[1]);
  Rescorrel = malloc(sizeof(float)*NumResidues*NumResidues);
  bzero(Rescorrel, sizeof(float)*NumResidues*NumResidues);
  
  // If rank 0, we'll load the data, distribute it to the other nodes, and write stuff to disk
  int ret;
  if(Rank == 0)
    ret = Master(argc, argv);
  else
    ret = Worker(argc, argv);
  
  MPI_Finalize();
  return ret;
}