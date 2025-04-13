#!/usr/bin/env bash

# Usage: ./install_polymv.sh <env_name> <installation_folder>
ENV_NAME=$1
INSTALLATION_FOLDER=$2

if [ -z "$ENV_NAME" ] || [ -z "$INSTALLATION_FOLDER" ]; then
  echo "❌ Usage: $0 <env_name> <installation_folder>"
  exit 1
fi

echo "🔧 Installing polymv..."

# Setup Conda
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Install polymv
pip install -e .

echo "✅ polymv installed."
