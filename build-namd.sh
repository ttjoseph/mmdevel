TOP=$PWD
if [ ! -d charm-6.8.2 ]; then
    tar xf charm-6.8.2.tar
fi

# Build and test the Charm++/Converse library (MPI version):
cd charm-6.8.2
env MPICXX=mpicxx ./build charm++ mpi-linux-x86_64 --with-production
# cd mpi-linux-x86_64/tests/charm++/megatest
# make pgm
# mpiexec -n 4 ./pgm   (run as any other MPI program on your cluster)
# cd ../../../../..

# Download and install TCL and FFTW libraries:
# (the precompiled FFTW provided by NAMD dev team is old and doesn't work for me)
cd $TOP
if [ ! -d fftw-3.3.8 ]; then
    wget http://www.fftw.org/fftw-3.3.8.tar.gz
    tar xf fftw-3.3.8.tar.gz
fi
mkdir -p fftw3
F=$PWD/fftw3
cd fftw-3.3.8
CFLAGS=-fPIC ./configure --enable-shared --enable-float --prefix=$F
make -j8 install
cd $TOP

# Download TCL libraries if necessary, in a janky fashion
if [ ! -d tcl ]; then
    wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.5.9-linux-x86_64.tar.gz
    tar xzf tcl8.5.9-linux-x86_64.tar.gz
    mv tcl8.5.9-linux-x86_64 tcl
fi
if [ ! -d tcl-threaded ]; then
    wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.5.9-linux-x86_64-threaded.tar.gz
    tar xzf tcl8.5.9-linux-x86_64-threaded.tar.gz
    mv tcl8.5.9-linux-x86_64-threaded tcl-threaded
fi

# Set up build directory and compile:
./config Linux-x86_64-g++ --with-fftw3 --fftw-prefix $PWD/fftw3 --charm-arch mpi-linux-x86_64
cd Linux-x86_64-g++
make -j8

# Quick tests using one and two processes (ethernet version):
#   (this is a 66-atom simulation so don't expect any speedup)
#   ./namd2
#   ./namd2 src/alanin
#   ./charmrun ++local +p2 ./namd2
#   ./charmrun ++local +p2 ./namd2 src/alanin
#   (for MPI version, run namd2 binary as any other MPI executable)

# Longer test using four processes:
#   wget http://www.ks.uiuc.edu/Research/namd/utilities/apoa1.tar.gz
#   tar xzf apoa1.tar.gz
#   ./charmrun ++local +p4 ./namd2 apoa1/apoa1.namd
#   (FFT optimization will take a several seconds during the first run.)
# Edit Make.charm to point at .rootdir/charm-6.8.2 or the full path to the
