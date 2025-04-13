#!/usr/bin/env bash

# Usage: ./install_dependencies.sh <pkg_manager> <env_name>
PKG_MANAGER=$1
ENV_NAME=$2

if [ -z "$PKG_MANAGER" ] || [ -z "$ENV_NAME" ]; then
  echo "❌ Usage: $0 <conda|mamba> <env_name>"
  exit 1
fi

if [[ "$PKG_MANAGER" != "conda" && "$PKG_MANAGER" != "mamba" ]]; then
  echo "❌ Invalid package manager: $PKG_MANAGER. Must be 'conda' or 'mamba'."
  exit 1
fi

echo "📦 Installing dependencies with ${PKG_MANAGER} in environment '${ENV_NAME}'..."

# Activate environment
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Core dependencies
${PKG_MANAGER} install \
  autoconf automake bison cmake cython doxygen flex gcc gfortran gxx help2man \
  libcurl libgcrypt libgfortran5 libtool m4 make pkg-config texinfo zlib -y

# Extra channel dependency (adjust as needed)
${PKG_MANAGER} install trung::libcheck -y

echo "✅ Dependencies installed using ${PKG_MANAGER}."
