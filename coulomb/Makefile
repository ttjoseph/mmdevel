CC=mpicc
CFLAGS=-O2
LDFLAGS=-lz -lm

all : coulomb correl

coulomb : coulomb.o 

correl : correl.o

%.o : %.c
	$(CC) $(CFLAGS) -c $<
	
clean :
	-rm -f coulomb coulomb.o
