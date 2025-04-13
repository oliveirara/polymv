	#!/usr/bin/env bash

# Usage: ./install_mpsolve.sh <pkg_manager> <env_name> <installation_folder>
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

echo "🔧 Installing MPSolve..."

# Setup Conda
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Install MPSolve from source
git clone https://github.com/robol/MPSolve.git
cd MPSolve

# Configure and install MPSolve to the specified directory
find . -type f -exec dos2unix {} +
bash autogen.sh
./configure --prefix=${INSTALLATION_FOLDER} CFLAGS='-Wall -march=native -O3' LDFLAGS="-L${INSTALLATION_FOLDER}/lib -L${CONDA_PREFIX}/lib" CPPFLAGS="-I${INSTALLATION_FOLDER}/include -I${CONDA_PREFIX}/include" --disable-examples --disable-documentation
find ./src/mpsolve -name 'Makefile' -exec sed -i 's/-lgmp -lm/-lgmp -lm -lgmpxx/g' {} +
make -j "$(nproc)"
make install
cd ..
rm -rf MPSolve

echo "✅ MPSolve installed."
