include $(GOROOT)/src/Make.$(GOARCH)

TARG=mpi

CGOFILES=\
	mpi.go 

CGO_CFLAGS=-I/usr/include/openmpi/1.2.3-gcc -m64 -pthread 
CGO_LDFLAGS=-L/usr/lib64/openmpi/1.2.3-gcc -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl

include $(GOROOT)/src/Make.pkg
