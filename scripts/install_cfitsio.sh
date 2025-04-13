	#!/usr/bin/env bash

# Usage: ./install_cfitsio.sh <pkg_manager> <env_name> <installation_folder>
PKG_MANAGER=$1
ENV_NAME=$2
INSTALLATION_FOLDER=$3

if [ -z "$PKG_MANAGER" ] || [ -z "$ENV_NAME" ] || [ -z "$INSTALLATION_FOLDER" ]; then
  echo "❌ Usage: $0 <conda|mamba> <env_name> <installation_folder>"
  exit 1
fi

if [[ "$PKG_MANAGER" != "conda" && "$PKG_MANAGER" != "mamba" ]]; then
  echo "❌ Invalid package manager: $PKG_MANAGER. Must be 'conda' or 'mamba'."
  exit 1
fi

echo "🔧 Installing CFITSIO..."

# Setup Conda
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Install CFITSIO from source
wget https://github.com/HEASARC/cfitsio/archive/refs/tags/cfitsio-4.5.0.zip
unzip cfitsio-4.5.0.zip
cd cfitsio-cfitsio-4.5.0

# Configure and install CFITSIO to the specified directory
./configure --prefix=${INSTALLATION_FOLDER} --enable-shared --disable-curl LDFLAGS="-L${INSTALLATION_FOLDER}/lib -L${CONDA_PREFIX}/lib" CPPFLAGS="-I${INSTALLATION_FOLDER}/include -I${CONDA_PREFIX}/include"
make -j "$(nproc)"
make install
cd ..
rm -rf cfitsio*

echo "✅ CFITSIO installed."
