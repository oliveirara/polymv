	#!/usr/bin/env bash

# Usage: ./install_chealpix.sh <env_name> <installation_folder>
ENV_NAME=$1
INSTALLATION_FOLDER=$2

if [ -z "$ENV_NAME" ] || [ -z "$INSTALLATION_FOLDER" ]; then
  echo "❌ Usage: $0 <env_name> <installation_folder>"
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
make shared CFITSIO_INCDIR=${CONDA_PREFIX}/include CFITSIO_LIBDIR=${CONDA_PREFIX}/lib
make install LIBDIR=${INSTALLATION_FOLDER}/lib INCDIR=${INSTALLATION_FOLDER}/include RANLIB="ar -ts"
cd ..
rm -rf chealpix

echo "✅ chealpix installed."
