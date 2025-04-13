	#!/usr/bin/env bash

# Usage: ./install_nlopt.sh <env_name> <installation_folder>
ENV_NAME=$1
INSTALLATION_FOLDER=$2

if [ -z "$ENV_NAME" ] || [ -z "$INSTALLATION_FOLDER" ]; then
  echo "❌ Usage: $0 <env_name> <installation_folder>"
  exit 1
fi

echo "🔧 Installing NLOPT..."

# Setup Conda
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Install NLOPT from source
wget https://github.com/stevengj/nlopt/archive/refs/tags/v2.9.0.tar.gz
tar -xvf v2.9.0.tar.gz
cd nlopt-2.9.0

# Configure and install NLOPT to the specified directory
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALLATION_FOLDER} ..
make -j "${nproc}"
make install
cd ../..
rm -rf v2.9.0.tar.gz nlopt-2.9.0

echo "✅ NLOPT installed."
