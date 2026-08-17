#!/usr/bin/env python3
"""
Figure 1 Generation Script: Transcriptomic Overview
Project: RNA-Seq Analysis of Salmonella enterica SL1344 (InSPI2 vs LSP)
Pipeline: READemption v2.0.4 + segemehl v0.3.4 + DESeq2 v1.50.2
"""

import os
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Paths to intermediate PDF-extracted panels
panel_dir = "/tmp/pdf_pngs"
img_pca = Image.open(os.path.join(panel_dir, "pca_heatmap-1.png"))
img_heatmap = Image.open(os.path.join(panel_dir, "pca_heatmap-2.png"))
img_volcano = Image.open(os.path.join(panel_dir, "volcano-1.png"))
img_ma = Image.open(os.path.join(panel_dir, "ma-1.png"))
img_scatter = Image.open(os.path.join(panel_dir, "scatter-02.png"))
img_align = Image.open(os.path.join(panel_dir, "align_reads-1.png"))

def autocrop_image(img, border=15):
    diff = ImageOps.invert(img.convert('RGB'))
    bbox = diff.getbbox()
    if bbox:
        left = max(0, bbox[0] - border)
        top = max(0, bbox[1] - border)
        right = min(img.width, bbox[2] + border)
        bottom = min(img.height, bbox[3] + border)
        return img.crop((left, top, right, bottom))
    return img

cropped_pca = autocrop_image(img_pca)
cropped_heatmap = autocrop_image(img_heatmap)
cropped_volcano = autocrop_image(img_volcano)
cropped_ma = autocrop_image(img_ma)
cropped_scatter = autocrop_image(img_scatter)
cropped_align = autocrop_image(img_align)

# Setup figure canvas (16:10 aspect ratio)
plt.rcParams['font.family'] = 'DejaVu Sans'
fig = plt.figure(figsize=(22, 14), dpi=300, facecolor='#ffffff')

# Main title banner
fig.text(0.5, 0.965, "RNA-Seq Transcriptomic Analysis of Salmonella enterica (SL1344)",
         ha='center', va='center', fontsize=23, fontweight='bold', color='#0f172a')
fig.text(0.5, 0.938, "SPI-2 Virulence Inducing (InSPI2) vs Control (LSP) | READemption & DESeq2 Pipeline",
         ha='center', va='center', fontsize=13.5, fontweight='medium', color='#334155')

# Define grid layout (2 rows x 3 columns)
gs = gridspec.GridSpec(2, 3, figure=fig,
                       left=0.04, right=0.96,
                       top=0.89, bottom=0.06,
                       wspace=0.18, hspace=0.26)

panels = [
    (gs[0, 0], cropped_pca, "A. Principal Component Analysis (PCA)", "Condition-specific sample separation (PC1: 78% var)"),
    (gs[0, 1], cropped_heatmap, "B. Sample Distance Matrix", "Hierarchical clustering of Euclidean distances"),
    (gs[0, 2], cropped_align, "C. Read Mapping Statistics", "Segemehl alignment distribution (Chromosome & Plasmids)"),
    (gs[1, 0], cropped_volcano, "D. Volcano Plot (InSPI2 vs LSP)", "Differentially expressed genes (-log10 padj vs log2FC)"),
    (gs[1, 1], cropped_ma, "E. MA Plot", "Log2 fold change vs mean normalized counts"),
    (gs[1, 2], cropped_scatter, "F. Expression Correlation", "Replicate expression concordance & dynamic range")
]

for slot, img, title, subtitle in panels:
    ax = fig.add_subplot(slot)
    ax.imshow(img)
    ax.axis('off')
    ax.set_title(f"{title}\n{subtitle}", fontsize=11.5, fontweight='bold', color='#0f172a', pad=10)

# Footer banner
footer_text = "Analysis Pipeline: READemption v2.0.4 • segemehl v0.3.4 • DESeq2 v1.50.2 • Organism: Salmonella enterica serovar Typhimurium SL1344"
fig.text(0.5, 0.022, footer_text, ha='center', va='center', fontsize=11, color='#64748b', fontweight='semibold')

# Output destinations
targets = [
    ("figures/Figure1_Transcriptomic_Overview.png", "figures/Figure1_Transcriptomic_Overview.pdf"),
    ("READemption_analysis/output/Figure1_Transcriptomic_Overview.png", "READemption_analysis/output/Figure1_Transcriptomic_Overview.pdf")
]

for png_p, pdf_p in targets:
    os.makedirs(os.path.dirname(png_p), exist_ok=True)
    plt.savefig(png_p, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    plt.savefig(pdf_p, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')

plt.close()
print("Figure 1 (Transcriptomic Overview) successfully generated.")
