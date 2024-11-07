# polyMV Compilation Guide 🚀

To compile the C code for polyMV, you must first install the following software:

## MPSolve 📚
[See here](../docs/MPSolve/README.md)

## Automake 1.17 🔧
```bash
wget https://ftp.gnu.org/gnu/automake/automake-1.17.tar.gz
tar -xvf automake-1.17.tar.gz 
cd automake-1.17/
./configure
make shared
make -j 10
sudo make install
sudo ldconfig
```

## CFITSIO 🌌
```bash
git clone https://github.com/HEASARC/cfitsio.git
cd cfitsio/
./configure
make shared
sudo make install
```

## CHealpix 🌍
```bash
git clone https://github.com/fabienbaron/chealpix.git
cd chealpix/
make CFITSIO_INCDIR=/usr/local/include CFITSIO_LIBDIR=/usr/local/lib
sudo make install LIBDIR=/usr/local/lib INCDIR=/usr/local/include RANLIB="ar -rsv"
sudo ldconfig
```

## nlopt 📈
```bash
wget https://github.com/stevengj/nlopt/archive/v2.7.1.tar.gz
tar -xvf v2.7.1.tar.gz 
cd nlopt-2.7.1/
cmake . && make && sudo make install
```

## Compilation Command 🛠️
After installing the required software, navigate to the src folder and run the following command to compile polyMV:

```bash
gcc polymv.c -o polyMV -Wall -march=native -O3 -fopenmp -lcfitsio -lmps -lm -lgmp -lchealpix -lcfitsio -lstdc++ -ffast-math -lnlopt
```

And that's it! 🎉
