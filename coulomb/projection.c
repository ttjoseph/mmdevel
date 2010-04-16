// C/pthreads version of the project-matrix function in correl.go
// Tom Joseph <thomas.joseph@mssm.edu>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#define NUM_THREADS 16
#define LINE_BUF_SIZE (2*1024*1024)

int NumPairs, NumResidues, *Pairs;
float *Correl, *Rescorrel;

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

// Does 1/NUM_THREADS part of the work of projecting the correlation matrix back onto a
// residue-residue matrix.
void* projectPiece(void *ptr) {
  int *data = (int*)ptr;
  int rank = data[0];
  int numNodes = data[1];
  printf("[%d] I am a thread. One of %d threads. But I am special.\n", rank, numNodes);
  
  int res_i = 0, res_j = 0;
  for(res_i = rank; res_i < NumResidues; res_i += numNodes) {
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
    float percent = (float)res_i;
    percent *= percent;
    percent /= (float) (NumResidues*NumResidues);
    printf("[%d] Finished row %d of %d, about %.0f%% done.\n", rank, res_i+1, NumResidues, percent*100);
  } // res_i
  
  
  return NULL;
}

int main (int argc, char const *argv[]) {
  printf("projection (C version) - T. Joseph <thomas.joseph@mssm.edu>\n\n");

  if(argc < 2) {
    printf("Usage: <num-residues>\n\nWe assume you want to load correl.txt and pairs.txt and write to rescorrel.txt.\n");
    return 1;
  }

  // Open and read the first line of correl.txt to decide on its dimensions
  FILE *fp = fopen("correl.txt", "r");
  char *buf = malloc(LINE_BUF_SIZE);
  fgets(buf, LINE_BUF_SIZE, fp);
  
  NumPairs = 1;
  strtok(buf, " "); // Assume there's at least one token
  while(strtok(NULL, " "))
    NumPairs++;
  
  NumResidues = atoi(argv[1]);
    
  printf("I think there are %d pairs represented in correl.txt, for a structure with %d residues.\n", NumPairs,
    NumResidues);
  rewind(fp);
  
  // Now we know how big the matrix is, so we can inhale it and begin processing
  printf("Allocating %.0fMB for the correl matrix. Loading...\n", ((double) NumPairs)*NumPairs*sizeof(float)/1048576);
  Correl = malloc(sizeof(float)*NumPairs*NumPairs);
  Pairs = malloc(sizeof(int)*NumPairs*2);
  // TODO: load Pairs
  Rescorrel = malloc(sizeof(float)*NumResidues*NumResidues);
  
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
    printf("ERROR: We didn't load as many numbers from correl.txt as we should have.\n");
    return 1;
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
  
  pthread_t threads[NUM_THREADS];
  int data[NUM_THREADS*2];
  
  // For N threads, the ith thread will get every Nth residue starting from residue i.
  int i;
  for(i = 0; i < NUM_THREADS; i++) {
    data[i*2] = i;
    data[i*2+1] = NUM_THREADS;
    pthread_create(&threads[i], NULL, projectPiece, &data[i*2]);
  }

  for(i = 0; i < NUM_THREADS; i++) {
    pthread_join(threads[i], NULL);
  }

  dumpFloatMatrixToFile("rescorrel.txt", Rescorrel, NumResidues, NumResidues);
  printf("Saved residue correlation matrix to rescorrel.txt.\n");
  
  return 0;
}
