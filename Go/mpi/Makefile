include $(GOROOT)/src/Make.$(GOARCH)

TARG=mpi

CGOFILES=\
	mpi.go 

# CGO_CFLAGS=-I/usr/include/openmpi/1.2.3-gcc -m64 -pthread 
# CGO_LDFLAGS=-L/usr/lib64/openmpi/1.2.3-gcc -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl

CGO_CFLAGS=-I/usr/mpi/intel/mvapich-1.1.0/include -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -D__INTEL_COMPILER -DMALLOC_RET_VOID=1 -Wl,-rpath-link -Wl,/usr/mpi/intel/mvapich-1.1.0/lib/shared
CGO_LDFLAGS=-L/usr/mpi/intel/mvapich-1.1.0/lib/shared -L/usr/mpi/intel/mvapich-1.1.0/lib -lmpich -L/usr/lib64 -libverbs -libumad -libcommon -lpthread -lpthread -lrt

include $(GOROOT)/src/Make.pkg
