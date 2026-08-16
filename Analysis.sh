#!/bin/bash
# ==============================================================================
# RNA-Seq Analysis Pipeline (READemption)
# Species: Salmonella enterica serovar Typhimurium
# Conditions: InSPI2 vs LSP
# ==============================================================================

set -e

# Set environment PATH
export PATH="$HOME/miniforge3/envs/reademption/bin:$HOME/miniforge3/bin:$PATH"

PROJECT_DIR="READemption_analysis"
PROCESSES=4

echo "=========================================================="
echo "Starting RNA-seq Analysis Pipeline with READemption"
echo "=========================================================="

# 1. Align Reads
echo "[Step 1/7] Aligning reads using segemehl..."
reademption align \
    --project_path "$PROJECT_DIR" \
    --processes "$PROCESSES" \
    --progress \
    --check_for_existing_files

# 2. Generate Coverage Files (Wiggle)
echo "[Step 2/7] Generating coverage (wiggle) files..."
reademption coverage \
    --project_path "$PROJECT_DIR" \
    --processes "$PROCESSES" \
    --check_for_existing_files

# 3. Gene-wise Quantification
echo "[Step 3/7] Performing gene-wise quantification..."
reademption gene_quanti \
    --project_path "$PROJECT_DIR" \
    --processes "$PROCESSES" \
    --check_for_existing_files

# 4. Differential Gene Expression Analysis with DESeq2
echo "[Step 4/7] Running differential expression analysis (InSPI2 vs LSP)..."
reademption deseq \
    --project_path "$PROJECT_DIR" \
    --libs "InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2" \
    --conditions "InSPI2,InSPI2,LSP,LSP" \
    --replicates "1,2,1,2" \
    --libs_by_species salmonella="InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2"

# 5. Visualizations - Alignments
echo "[Step 5/7] Generating alignment visualisations..."
reademption viz_align \
    --project_path "$PROJECT_DIR"

# 6. Visualizations - Gene Quantification
echo "[Step 6/7] Generating gene quantification visualisations (PCA, Heatmaps)..."
reademption viz_gene_quanti \
    --project_path "$PROJECT_DIR"

# 7. Visualizations - DESeq2 Differential Expression
echo "[Step 7/7] Generating DESeq2 visualisations (Volcano plots, MA plots)..."
reademption viz_deseq \
    --project_path "$PROJECT_DIR"

echo "=========================================================="
echo "RNA-seq Analysis Pipeline successfully completed!"
echo "=========================================================="
