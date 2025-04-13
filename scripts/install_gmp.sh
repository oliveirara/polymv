#!/usr/bin/env bash

# Usage: ./install_gmp.sh <env_name> <installation_folder>
ENV_NAME=$1
INSTALLATION_FOLDER=$2

if [ -z "$ENV_NAME" ] || [ -z "$INSTALLATION_FOLDER" ]; then
  echo "❌ Usage: $0 <env_name> <installation_folder>"
  exit 1
fi

echo "🔧 Installing GMP..."

# Setup Conda
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Install GMP from source
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar -xvf gmp-6.3.0.tar.xz
cd gmp-6.3.0

# Configure and install GMP to the specified directory
./configure --prefix="${INSTALLATION_FOLDER}" --enable-cxx
make -j "$(nproc)"
make install
cd ..
rm -rf gmp-6.3.0*

echo "✅ GMP installed."
