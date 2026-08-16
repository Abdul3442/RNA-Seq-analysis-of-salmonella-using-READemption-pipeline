#!/bin/bash
# ==============================================================================
# RNA-Seq Analysis Environment Setup (READemption Pipeline)
# ==============================================================================

set -e

# Step 1: Install Miniforge if conda/mamba is not available
if ! command -v mamba &> /dev/null && ! command -v conda &> /dev/null; then
    echo "Installing Miniforge..."
    curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
    rm /tmp/miniforge.sh
    export PATH="$HOME/miniforge3/bin:$PATH"
fi

# Step 2: Configure conda channels
conda config --add channels conda-forge
conda config --add channels bioconda

# Step 3: Create the 'reademption' environment with Python 3.9 and core tools
mamba create -n reademption python=3.9 -y
mamba install -n reademption -c bioconda -c conda-forge segemehl samtools -y
mamba install -n reademption -c conda-forge -c bioconda bioconductor-deseq2 r-base -y

# Step 4: Install READemption Python package
"$HOME/miniforge3/envs/reademption/bin/pip" install READemption

echo "READemption environment successfully set up!"
