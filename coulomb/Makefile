MPICC=mpicc
CC=icc
CFLAGS=-O2 -g -Wall -std=c99
LDFLAGS=-lz -lm -lpthread

all : coulomb correl projection

coulomb : coulomb.o 
	$(MPICC) coulomb.o $(LDFLAGS) -o coulomb

correl : correl.o
	$(MPICC) correl.o $(LDFLAGS) -o correl

projection : projection.c
	$(CC) projection.c $(CFLAGS) $(LDFLAGS) -o projection

%.o : %.c
	$(MPICC) $(CFLAGS) -c $<
	
.PHONY : clean

clean :
	-rm -f coulomb coulomb.o correl correl.o projection projection.o
