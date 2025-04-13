#!/usr/bin/env bash

# Usage: ./create_env.sh <pkg_manager> <env_name> <python_version>
PKG_MANAGER=$1
ENV_NAME=$2
PYTHON_VERSION=$3

if [ -z "$PKG_MANAGER" ] || [ -z "$ENV_NAME" ]; then
  echo "❌ Usage: $0 <conda|mamba> <env_name>"
  exit 1
fi

echo "🐍 Creating environment '${ENV_NAME}' with Python ${PYTHON_VERSION} using ${PKG_MANAGER}..."

# Load Conda shell functions
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Create the environment
${PKG_MANAGER} create -n "${ENV_NAME}" python=${PYTHON_VERSION} -y

echo "✅ Environment '${ENV_NAME}' created with ${PKG_MANAGER}."