MPICC=mpicc
CC=icc
CFLAGS=-O2 -g -Wall -std=c99
LDFLAGS=-lz -lm -lpthread

all : coulomb correl projection projection-mpi

coulomb : coulomb.o 
	$(MPICC) coulomb.o $(LDFLAGS) -o coulomb

correl : correl.o
	$(MPICC) correl.o $(LDFLAGS) -o correl

projection-mpi : projection-mpi.o
	$(MPICC) projection-mpi.o $(LDFLAGS) -o projection-mpi

projection : projection.c
	$(CC) projection.c $(CFLAGS) $(LDFLAGS) -o projection

%.o : %.c
	$(MPICC) $(CFLAGS) -c $<
	
.PHONY : clean

clean :
	-rm -f coulomb coulomb.o correl correl.o projection projection.o
