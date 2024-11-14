# polyMV Compilation Guide

Follow these steps to compile the polyMV program. If you have `sudo` permissions, follow `Installation 1`. Otherwise, follow `Installation 2`. This README where prepared using a Docker image from Ubuntu distribution.

## Installation 1

### System packages

Update and install the packages below:

```bash
sudo apt update
sudo apt install autoconf automake autotools-dev bison check cmake curl flex g++ gcc gfortran git help2man libgmp-dev libmps-dev libpthread-stubs0-dev libtool m4 make pkg-config texinfo wget
```

### MPSolve

Please refer to the [official installation steps](https://numpi.dm.unipi.it/scientific-computing-libraries/mpsolve/) to install MPSolve. However, in some cases, you may need to build MPSolve on your machine, particularly if you're using a cluster where you don't have `sudo` permissions. The steps below will guide you through building MPSolve from source. If you have any questions, feel free to let us know!

```bash
git clone https://github.com/robol/MPSolve.git
cd MPSolve/
bash autogen.sh
./configure
make -j $(nproc)
make check
sudo make install
sudo ldconfig
```

Verify the Installation. Reset the shell and check the MPSolve version:

```bash
mpsolve -v
```

### Automake

```bash
wget https://ftp.gnu.org/gnu/automake/automake-1.17.tar.gz
tar -xvf automake-1.17.tar.gz
cd automake-1.17/
./configure
make shared
make -j $(nproc)
sudo make install
sudo ldconfig
```

### CFITSIO

```bash
git clone https://github.com/HEASARC/cfitsio.git
cd cfitsio/
./configure --enable-shared
make -j $(nproc)
sudo make install
sudo ldconfig
```

### CHealpix

```bash
git clone https://github.com/fabienbaron/chealpix.git
cd chealpix/
make shared CFITSIO_INCDIR=/usr/local/include CFITSIO_LIBDIR=/usr/local/lib
sudo make install LIBDIR=/usr/local/lib INCDIR=/usr/local/include RANLIB="ar -ts"
sudo ldconfig
```

### HDF5

```bash
wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
tar -xvf hdf5-1.14.5.tar.gz
cd hdf5-1.14.5
./configure --enable-fortran --enable-cxx 
make -j $(nproc)
sudo make install
sudo ldconfig
```

### NLOPT

```bash
wget https://github.com/stevengj/nlopt/archive/refs/tags/v2.9.0.tar.gz
tar -xvf v2.9.0.tar.gz
cd nlopt-2.9.0
mkdir build
cd build
cmake ..
make -j $(nproc)
sudo make install
sudo ldconfig
```

### polyMV

After installing the required software:

```bash
git clone https://github.com/oliveirara/polyMV.git -b polyMV_c
cd polyMV/src
gcc polymv.c -o polyMV -Wall -march=native -O3 -fopenmp -lcfitsio -lgmp -lgmpxx -lm -lmps -lchealpix -lstdc++ -ffast-math -lnlopt -lhdf5
```

And that's it! 🎉

## Installation 2

Sometimes you don't have `sudo` powers. However, if you specify correctly the PATH, you compile and install whatever library/software. First, let's define where we are going to install most of all necessary packages:

```bash
mkdir ~/software
export INSTALLATION_FOLDER=/home/${USER}/software/
```

### Install Miniforge

Download and install Miniforge:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh 
source .bashrc
```

### Install packages using mamba

Optional: you can create your own environment to install all packages. For this tutorial, we will install all packages on `(base)`:

```bash
mamba install gcc gxx autoconf cmake gmp automake bison flex libtool m4 cython help2man libgfortran5 pkg-config texinfo doxygen make gfortran zlib libgcrypt libcurl zlib -y
mamba install trung::libcheck -y
source .bashrc
```

### GMP

```bash
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar -xvf gmp-6.3.0.tar.xz 
cd gmp-6.3.0
./configure --prefix=${INSTALLATION_FOLDER}/gmp --enable-cxx
make -j $(nproc)
make install
cd
```

### MPSolve

```bash
git clone https://github.com/robol/MPSolve.git
cd MPSolve/
bash autogen.sh
./configure --prefix=${INSTALLATION_FOLDER}/mpsolve CFLAGS='-Wall -march=native -O3' LDFLAGS="-L${INSTALLATION_FOLDER}/gmp/lib" CPPFLAGS="-I${INSTALLATION_FOLDER}/gmp/include" --disable-examples
find ./src/mpsolve -name 'Makefile' -exec sed -i 's/-lgmp -lm/-lgmp -lm -lgmpxx/g' {} +
make -j $(nproc)
make install
cd
```

### CFITSIO

```bash
git clone https://github.com/HEASARC/cfitsio.git
cd cfitsio/
./configure --prefix=${INSTALLATION_FOLDER}/cfitsio --enable-shared --disable-curl LDFLAGS="-L/home/${USER}/miniforge3/lib" CPPFLAGS="-I/home/${USER}/miniforge3/include"
make -j $(nproc)
make install
cd
```

If you create a new environment, the LDFLAGS and CPPFLAGS should be something like this: `"/home/${USER}/miniforge3/envs/${ENV_NAME}/"`.

### Chealpix

```bash
cd chealpix
make shared CFITSIO_INCDIR=${INSTALLATION_FOLDER}/cfitsio/include CFITSIO_LIBDIR=${INSTALLATION_FOLDER}/cfitsio/lib
mkdir ${INSTALLATION_FOLDER}/chealpix ${INSTALLATION_FOLDER}/chealpix/lib ${INSTALLATION_FOLDER}/chealpix/include
make install LIBDIR=${INSTALLATION_FOLDER}/chealpix/lib INCDIR=${INSTALLATION_FOLDER}/chealpix/include RANLIB="ar -ts"
cd
```

### HDF5

```bash
wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
tar -xvf hdf5-1.14.5.tar.gz
cd hdf5-1.14.5
./configure --prefix=${INSTALLATION_FOLDER}/hdf5 --enable-fortran --enable-cxx 
make -j $(nproc)
make install
cd
```

### NLOPT

```bash
wget https://github.com/stevengj/nlopt/archive/refs/tags/v2.9.0.tar.gz
tar -xvf v2.9.0.tar.gz
cd nlopt-2.9.0
mkdir build
cd build
mkdir ${INSTALLATION_FOLDER}/nlopt
cmake -DCMAKE_INSTALL_PREFIX=${INSTALLATION_FOLDER}/nlopt ..
make -j $(nproc)
make install
cd
```

### polyMV

Export some paths:

```bash
```

It will be convenient to export this file to `~/.bashrc` file:

```
export LD_LIBRARY_PATH=${INSTALLATION_FOLDER}/mpsolve/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${INSTALLATION_FOLDER}/chealpix/lib:$LD_LIBRARY_PATH
```

And compile polyMV:

 ```bash
git clone https://github.com/oliveirara/polyMV.git -b polyMV_c
cd polyMV/src
gcc polymv.c -o polyMV -Wall -march=native -O3 -fopenmp -I${INSTALLATION_FOLDER}/cfitsio/include -L${INSTALLATION_FOLDER}/cfitsio/lib -lcfitsio -I${INSTALLATION_FOLDER}/gpm/include -L${INSTALLATION_FOLDER}/gpm/lib -lgmp -lgmpxx -I${INSTALLATION_FOLDER}/mpsolve/include -L${INSTALLATION_FOLDER}/mpsolve/lib -lmps -lm -I${INSTALLATION_FOLDER}/chealpix/include -L${INSTALLATION_FOLDER}/chealpix/lib -lchealpix -lstdc++ -ffast-math -I${INSTALLATION_FOLDER}/nlopt/include -L${INSTALLATION_FOLDER}/nlopt/lib -lnlopt -I${INSTALLATION_FOLDER}/hdf5/include -L${INSTALLATION_FOLDER}/hdf5/lib -lhdf5
```

And that's it! 🎉
