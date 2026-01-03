
# Integrated proteomic and single-cell transcriptomic profiling elucidates immunomodulatory effects of L-serine in autism spectrum disorder
## 📌 Overview
This repository contains the computational workflow for analyzing single-cell RNA sequencing (scRNA-seq) data from PBMC samples derived from Autism Spectrum Disorder (ASD) patients treated with L-serine compared to Healthy Controls (HC).

The analysis pipeline focuses on:

- Global Profiling: Characterizing PBMC subsets and cell type composition changes.

- Gene Scoring: Evaluating ASD-risk gene expression (SFARI, TADA) and inflammatory signatures.

- CD4+ T Cell Dynamics: Trajectory inference (Slingshot), lineage-based differential expression (tradeSeq), and identification of "rescue" genes reversed by treatment.

- Network Analysis: Co-expression network analysis (hdWGCNA) on CD8+ T cells.

- Cell-Cell Communication: Interaction analysis using CellPhoneDB v5.
-------
##### Author:  *Jaemyung Jang*, Seungeun Yeo, Jong Pil Kim, Soonbong Baek, Hyun Jin Jung, Su-Kyeong Hwang, Youngshik Choe
-------

🛠 Prerequisites & Installation
### R
R Environment (v4.5.2)
Key R packages required for the analysis:

##### Core & Visualization 
```
install.packages(c("Seurat", "dplyr", "ggplot2", "ggpubr", "pheatmap", "viridis", "ggsci"))
```
##### Bioconductor
```
BiocManager::install(c("Slingshot", "tradeSeq", "UCell", "Nebulosa", "dittoSeq", "clusterProfiler", "ReactomePA", "org.Hs.eg.db"))
```
##### Networks
```
devtools::install_github("smorabit/hdWGCNA") 
library(hdWGCNA)
library(WGCNA)
```
### Python Environment (v3.9+)

Install the required Python packages using pip:
```
python -m pip install Plaintext pandas anndata  matplotlib cellphonedb>=5.0.0 ktplotspy
```

## 🚀 Workflow Instructions
### Step 1: Global PBMC Analysis (R)
 Run the initial profiling script to generate global UMAPs, density plots for markers, and calculate inflammatory scores.

+ Key Outputs:
- Figure 2B: UMAP by Group (Healthy vs ASD).
- Figure 2C: Density plots (Nebulosa) for lineage markers.
- Figure 2D: Cell type proportion bars.
- Figure 2E: SFARI/TADA gene set scoring.
- Supp Fig 2B: Inflammatory cytokine scoring.

### Step 2: CD4+ T Cell Trajectory & Rescue Analysis (R)
 Perform deep-dive analysis on CD4+ T cells, including trajectory inference and identifying genes "rescued" by L-serine treatment (Genes UP in ASD -> DOWN in Treat, and vice versa).

+ Key Outputs:
- Figure 3A: CD4+ T cell UMAP.
- Figure 3E-F: Rescue effect heatmaps and bar plots.
- Figure 4A: Slingshot trajectories (ASD vs Treat).
- Figure 4C: Smoothed gene expression heatmaps along pseudotime.
- Figure 4D: Individual gene trend plots.

### Step 3: Cell-Cell Interaction Analysis (Python)
 Execute the Python script to analyze interactions using CellPhoneDB v5.

1. Downloads CellPhoneDB v5 database.
2. Analyzes BO (Before OM / ASD) and AO (After OM / Treat) conditions.
3. Compares interaction counts between groups.

+ Key Outputs:
  - Figure5a1/2: Heatmaps of significant interaction counts.
  - Figure5b1/2: Dotplots of ligand-receptor pairs (CD4+ T focused).

### Step 4: CD4+ T Cell Network Analysis & Data Export (R)
Run hdWGCNA on CD8+ T cells to identify co-expression modules and export count matrices for CellPhoneDB.

+ Key Actions:
- Construct co-expression networks.
- Identify Differential Module Expression (DME).
- Export: Generates .txt count matrices and .tsv metadata files in results/v5.0.0/ for the Python step.


### ⚠️ Important Notes
- RDS Files: The input Seurat object (ASD_PBMC_subset_revision.RDS) is large and expected to be present in the data/ folder.

- WGCNA Threads: The WGCNA analysis requires significant RAM and CPU threads. Ensure allowWGCNAThreads() is configured for your machine.

### 📜 Citation
"under review"
