	#!/usr/bin/env bash

# Usage: ./install_chealpix.sh <pkg_manager> <env_name> <installation_folder>
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

echo "🔧 Installing chealpix..."

# Setup Conda
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Install chealpix from source
git clone https://github.com/fabienbaron/chealpix.git
cd chealpix

# Configure and install chealpix to the specified directory
make shared CFITSIO_INCDIR=${INSTALLATION_FOLDER}/include CFITSIO_LIBDIR=${INSTALLATION_FOLDER}/lib
make install LIBDIR=${INSTALLATION_FOLDER}/lib INCDIR=${INSTALLATION_FOLDER}/include RANLIB="ar -ts"
cd ..
rm -rf chealpix

echo "✅ chealpix installed."
