# MPSolve Installation Steps

Please refer to the [official installation steps](https://numpi.dm.unipi.it/scientific-computing-libraries/mpsolve/) to install MPSolve. However, in some cases, you may need to build MPSolve on your machine, particularly if you're using a cluster where you don't have `sudo` permissions. The steps below will guide you through building MPSolve from source. If you have any questions, feel free to let us know!

## 🚀 MPSolve Installation Guide

This readme where prepared using a Docker image from Ubuntu distribution. Follow these steps to install **MPSolve** assuming you are `root` or an user with `sudo` powers:

### 1. Install Miniforge 🛠️

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

### 2. (Optional) Create an Isolated Environment 🧑‍🔬

```bash
conda create -n mpsolve
```

### 3. Activate the Environment ⚙️

```bash
conda activate mpsolve
```

### 4. Install Required Packages 📦

If you're root:


```bash
# apt update
# apt install autoconf automake autotools-dev bison check cmake curl flex g++ gcc gfortran git help2man libgmp-dev libmps-dev libpthread-stubs0-dev libtool m4 make pkg-config texinfo wget
```

If you're user with `sudo` powers, just add `sudo` at the beggining of the line above (e.g., `sudo apt ...`).

### 5. Clone the MPSolve Repository 🧩

```bash
git clone https://github.com/robol/MPSolve.git
cd MPSolve/
```

### 6. Run the `autogen` Script 🔧

```bash
bash autogen.sh
```

### 7. Configure Installation Path 📁

Set the installation folder by running:

The default installation folder is `/usr/local/bin/mpsolve`. Then

```bash
./configure
```

Then compile and install:

```bash
make -j 10
make check
sudo make install
sudo ldconfig
```

### 8. Verify the Installation ✅

Reset the shell and check the MPSolve version:

```bash
mpsolve -v
```

## 🚀 Installation Guide for `che` Cluster Users

Follow these steps to install **MPSolve** on the `che` cluster:

### 1. Install Miniforge 🛠️

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

### 2. (Optional) Create an Isolated Environment 🧑‍🔬

```bash
conda create -n mpsolve
```

### 3. Activate the Environment ⚙️

```bash
conda activate mpsolve
```

### 4. Clone the MPSolve Repository 🧩

```bash
git clone https://github.com/robol/MPSolve.git
cd MPSolve/
```

### 5. Load Necessary Modules 🧰

To load the necessary modules, add the following lines to your `~/.bashrc` file or run these commands directly in your terminal:

```bash
module load autoconf/2.72-gcc-5.3.0
module load automake/1.16.5-gcc-5.3.0
module load autotools
module load bison/3.8.2-gcc-12.2.0
module load cmake/3.24.2-gcc-12.2.0
module load flex/2.6.3-gcc-12.2.0
module load gcc/5.3.0
module load gmp/6.2.1-gcc-5.3.0
module load libtool/2.4.7-gcc-5.3.0
module load m4/1.4.19-gcc-5.3.0
```

### 6. Install Required Packages 📦

```bash
conda install cython help2man libgfortran5 pkg-config texinfo
conda install trung::libcheck
```

### 7. Run the `autogen` Script 🔧

```bash
bash autogen.sh
```

### 8. Configure Installation Path 📁

Set the installation folder by running the `configure` script with the `--prefix` option. Replace `<your installation folder>` with your desired installation path, for example, `/home/<your-username>/software/mpsolve`:

```bash
./configure --prefix=/home/<your-username>/software/mpsolve
```

Then compile and install:

```bash
make -j 10
make check
make install
```

### 9. Verify the Installation ✅

Check the MPSolve version:

```bash
/home/<your-username>/software/mpsolve/bin/mpsolve -v
```

### 10. Add MPSolve to Your `PATH` 🛤️

Add the following line to your `~/.bashrc` file to ensure the MPSolve binary folder is in your `PATH`:

```bash
echo "export PATH=/home/<your-username>/software/mpsolve/bin:\${PATH}" >> ~/.bashrc
```

### 11. Reload Your Shell 🔄

Reload your shell, and you'll be able to run MPSolve from any directory:

```bash
source ~/.bashrc
```

🎉 **You're all set to use MPSolve!**
