CC=mpicc
CFLAGS=-O2
LDFLAGS=-lz

coulomb : coulomb.o 

%.o : %.c
	$(CC) $(CFLAGS) -c $<
	
clean :
	-rm -f coulomb coulomb.o
