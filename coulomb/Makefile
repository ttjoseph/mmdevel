CC=mpicc
CFLAGS=-O2 -g -Wall
LDFLAGS=-lz -lm

all : coulomb correl

coulomb : coulomb.o 

correl : correl.o

%.o : %.c
	$(CC) $(CFLAGS) -c $<
	
clean :
	-rm -f coulomb coulomb.o correl correl.o
