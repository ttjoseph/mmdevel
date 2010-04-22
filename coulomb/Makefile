CC=mpicc
CFLAGS=-O2 -g -Wall
LDFLAGS=-lz -lm -lpthread

all : coulomb correl projection

coulomb : coulomb.o 

correl : correl.o

projection : projection.o

%.o : %.c
	$(CC) $(CFLAGS) -c $<
	
clean :
	-rm -f coulomb coulomb.o correl correl.o projection projection.o
