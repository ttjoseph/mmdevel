CC=mpicc
CFLAGS=-O2 -g -Wall -std=c99
LDFLAGS=-lz -lm -lpthread

all : coulomb correl projection

coulomb : coulomb.o 
	$(CC) coulomb.o $(LDFLAGS) -o coulomb

correl : correl.o
	$(CC) correl.o $(LDFLAGS) -o correl

projection : projection.o
	$(CC) projection.o $(LDFLAGS) -o projection

%.o : %.c
	$(CC) $(CFLAGS) -c $<
	
.PHONY : clean

clean :
	-rm -f coulomb coulomb.o correl correl.o projection projection.o
