################################################################################
# Title: Integrated proteomic and single-cell transcriptomic profiling elucidates 
#        immunomodulatory effects of L-serine in autism spectrum disorder
# Subtitle: CD4+ T Cell Trajectory & Differential Expression Analysis
# Description: This script performs trajectory inference using Slingshot and lineage-based 
#   differential expression analysis using tradeSeq on CD4+ T cell subsets.
#   It visualizes gene expression trends along pseudotime for ASD and Treatment groups.
# Author: Jaemyung Jang / Korea Brain Research institue (piloter2@kbri.re.kr)
# Date: 2026-01-02
# R Version: 4.5.2
################################################################################

# ==============================================================================
# 1. Setup & Libraries
# ==============================================================================
suppressPackageStartupMessages({
  # Data Manipulation & Visualization
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(lattice)
  library(readxl)
  library(data.table)
  library(ggvenn)
  library(pheatmap)
  library(viridis)
  
  # Single-cell Analysis
  library(Seurat)
  library(Nebulosa)
  library(UCell)
  library(SingleCellExperiment)
  library(dittoSeq)
  
  # Trajectory & DE Analysis
  library(slingshot)
  library(tradeSeq)
})
# ==============================================================================
# 2. Data Loading & Preprocessing
# ==============================================================================
# --- Parameters & Paths ---
# NOTE: Users should adjust these paths to match their local environment.
BASE_DIR    <- getwd() 
DATA_DIR    <- file.path(BASE_DIR, "data")
RESULTS_DIR <- file.path(BASE_DIR, "results")
RDS_PATH    <- file.path(DATA_DIR, "ASD_PBMC_subset_revision.RDS")

# Load Seurat Object
if (!file.exists(RDS_PATH)) stop("RDS file not found. Check the path.")
PBMC.subset <- readRDS(RDS_PATH)

# Define Color Palette
custom_cols <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
curve_cols  <- c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "lightgrey", "darkgrey")

# ==============================================================================
# 3. Preprocessing & Subsetting
# ==============================================================================
# Initial SCT normalization (if needed before subsetting)
PBMC.subset <- PrepSCTFindMarkers(PBMC.subset, assay = "SCT")

# Subset for specific donors (Healthy Controls and specific Human samples)
target_donors <- c("hc_1", "hc_2", "hc_7", "hc_8", "human1", "human2")
PBMC.subsets <- subset(PBMC.subset, subset = Donor %in% target_donors)

# Re-run PrepSCT needed after subsetting for correct SCT model usage
PBMC.subsets <- PrepSCTFindMarkers(PBMC.subsets, assay = "SCT")

# ==============================================================================
# 4. CD4+ T Cell Subsetting and Visualization
# ==============================================================================

# Subset CD4+ T cells
CD4T <- subset(PBMC.subsets, subset = cell.type.coarse == "CD4+ T")

# Define Colors
cols_ <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
CD4T_sce <- as.SingleCellExperiment(CD4T, assay = "SCT")

# ==============================================================================
# 5. Trajectory Inference (Slingshot)
# ==============================================================================
# Run Slingshot
# Note: 'approx_points' reduces computational load for large datasets
sce_slingshot <- slingshot(
  CD4T_sce, 
  reducedDim = 'UMAP', 
  clusterLabels = colData(CD4T_sce)$cell.type.fine, 
  start.clus = "CD4n T-A", 
  approx_points = 250, 
  maxit = 1000, 
  thresh = 0.001
)

# Extract Lineages and Curves
lin <- getLineages(sce_slingshot, reducedDim = 'UMAP', clusterLabels = colData(CD4T_sce)$cell.type.fine, start.clus = "CD4n T-A")
crv <- getCurves(lin)

# Visualize Trajectories (ASD vs Treatment)
# Adding grouping info for visualization
colData(sce_slingshot)$cell.type <- colData(sce_slingshot)$cell.type.fine
colData(sce_slingshot)$id_new    <- colData(CD4T_sce)$id_new # Ensure metadata is consistent

# Figure 4a : ASD Condition
asd_sce <- sce_slingshot
colData(asd_sce)$cell.type[which(colData(asd_sce)$id_new != "2. Autism")] <- "BG"
p_traj_asd <- dittoDimPlot(
  asd_sce,
  var = "cell.type", size = 0.25, theme = theme_void(),
  reduction.use = "UMAP", trajectory.arrow.size = 0.2,
  add.trajectory.curves = slingCurves(sce_slingshot)[c(1, 2)],
  colors = c(32, 1, 6:7, 3, 5, 4, 2),
  legend.show = TRUE,
  main = paste0("ASD: ", sum(colData(sce_slingshot)$id_new == "2. Autism"), " cells"),
  sub = "Start cluster: CD4 T naive IL7R low"
)

# Figure 4a : Treatment (AST) Condition
ast_sce <- sce_slingshot
colData(ast_sce)$cell.type[which(colData(ast_sce)$id_new != "3. Treat")] <- "BG"
p_traj_ast <- dittoDimPlot(
  ast_sce,
  var = "cell.type", size = 0.25, theme = theme_void(),
  reduction.use = "UMAP", trajectory.arrow.size = 0.2,
  add.trajectory.curves = slingCurves(sce_slingshot)[c(1, 2)],
  colors = c(32, 1, 6:7, 3, 5, 4, 2),
  legend.show = TRUE,
  main = paste0("AST: ", sum(colData(sce_slingshot)$id_new == "3. Treat"), " cells"),
  sub = "Start cluster: CD4 T naive IL7R low"
)

# Save Trajectory Plots (Optional)
# ggsave(file.path(RESULTS_DIR, "Trajectory_ASD.pdf"), p_traj_asd)
# ggsave(file.path(RESULTS_DIR, "Trajectory_AST.pdf"), p_traj_ast)

# ==============================================================================
# 6. Differential Expression Analysis (tradeSeq)
# ==============================================================================

message("Running fitGAM (General Additive Models). This may take time...")
base::gc() # Garbage collection to free memory

# Determine optimal number of knots (nknots) based on previous analysis or AIC
sce_gam <- fitGAM(
  counts = counts(CD4T_sce), 
  sds = SlingshotDataSet(sce_slingshot), 
  conditions = factor(colData(CD4T_sce)$id_new),
  parallel = TRUE,
  BPPARAM = BiocParallel::bpparam(),
  nknots = 9, 
  verbose = TRUE
)

# Perform Association Test
rowData(sce_gam)$assocRes <- associationTest(sce_gam, lineages = TRUE, l2fc = log2(0.25))
assocRes <- rowData(sce_gam)$assocRes

# Rename columns for clarity
colnames(assocRes) <- gsub("2. Autism", "_ASD_before_OM", colnames(assocRes))
colnames(assocRes) <- gsub("3. Treat", "_ASD_after_OM", colnames(assocRes))

# Identify Significant Genes (FDR <= 0.01)
# Lineage 1
ASDGenes_L1 <- rownames(sce_gam)[which(p.adjust(assocRes$'pvalue_lineage1_condition_ASD_before_OM', "fdr") <= 0.01)]
ASTGenes_L1 <- rownames(sce_gam)[which(p.adjust(assocRes$'pvalue_lineage1_condition_ASD_after_OM', "fdr") <= 0.01)]

# Lineage 2
ASDGenes_L2 <- rownames(sce_gam)[which(p.adjust(assocRes$'pvalue_lineage2_condition_ASD_before_OM', "fdr") <= 0.01)]
ASTGenes_L2 <- rownames(sce_gam)[which(p.adjust(assocRes$'pvalue_lineage2_condition_ASD_after_OM', "fdr") <= 0.01)]

# ==============================================================================
# 7. Overlap with EV Genes & Visualization
# ==============================================================================

# Load EV-associated gene lists
ev_pos_path <- file.path(DATA_DIR, "ASD_PBMC_EV-pos_genes.txt")
ev_neg_path <- file.path(DATA_DIR, "ASD_PBMC_EV-neg_genes.txt")

if(file.exists(ev_pos_path) & file.exists(ev_neg_path)){
  EVgenes_pos <- fread(ev_pos_path)[,1] %>% pull()
  EVgenes_neg <- fread(ev_neg_path)[,1] %>% pull()
  EVgenes_all <- c(EVgenes_pos, EVgenes_neg)
  
  # Venn Diagrams: Overlap between trajectory-associated genes and EV genes
  # Lineage 1
  venn_l1 <- ggvenn(
    list("ASD (Before)" = ASDGenes_L1, "AST (After)" = ASTGenes_L1, "EV Genes" = EVgenes_all), 
    fill_color = c("red", "blue", "black")
  )
  # Figure 4b  
  print(venn_l1)
  
  # Lineage 2
  venn_l2 <- ggvenn(
    list("ASD (Before)" = ASDGenes_L2, "AST (After)" = ASTGenes_L2, "EV Genes" = EVgenes_all), 
    fill_color = c("red", "blue", "black")
  )
  # Figure 4b  
  print(venn_l2)
  
} else {
  warning("EV gene lists not found. Skipping Venn diagrams.")
  EVgenes_all <- character(0) # Empty vector to prevent errors downstream if files missing
}

# ==============================================================================
# 8. Heatmap Generation (Smoothed Expression)
# ==============================================================================

message("Generating smoothed expression heatmaps...")

# Define helper to calculate predicted smoothers for gene intersections
get_smoothers <- function(genes, model, n_points = 100) {
  if(length(genes) == 0) return(NULL)
  predictSmooth(model, gene = genes, nPoints = n_points, tidy = FALSE)
}

# Calculate intersections with EV genes
target_genes_L1_ASD_specific <- intersect(EVgenes_all, setdiff(ASDGenes_L1, ASTGenes_L1))
target_genes_L1_AST_specific <- intersect(EVgenes_all, setdiff(ASTGenes_L1, ASDGenes_L1))
target_genes_L1_Common       <- intersect(EVgenes_all, intersect(ASTGenes_L1, ASDGenes_L1))

# Predict Smoothers (Lineage 1)
yhat_L1_ASD_spec <- get_smoothers(target_genes_L1_ASD_specific, sce_gam)
yhat_L1_AST_spec <- get_smoothers(target_genes_L1_AST_specific, sce_gam)
yhat_L1_Common   <- get_smoothers(target_genes_L1_Common, sce_gam)

# Function to draw and save heatmaps
save_heatmap <- function(mat, filename, title = "") {
  if (is.null(mat)) return(NULL)
  scaled_mat <- t(scale(t(mat))) # Z-score scaling across pseudotime
  
  pheatmap(
    scaled_mat,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    border_color = "white",
    main = title,
    filename = file.path(RESULTS_DIR, filename)
  )
}

# Save Heatmaps for Lineage 1
# Figure 4c
save_heatmap(yhat_L1_AST_spec, "Figure4C_AST.pdf")
save_heatmap(yhat_L1_ASD_spec, "Figure4C_ASD.pdf")
save_heatmap(yhat_L1_Common,   "Figure4C_common.pdf")

# Combined Heatmap (Track A)
if (!is.null(yhat_L1_AST_spec) && !is.null(yhat_L1_ASD_spec) && !is.null(yhat_L1_Common)) {
  combined_mat <- rbind(
    t(scale(t(yhat_L1_AST_spec))),
    t(scale(t(yhat_L1_ASD_spec))),
    t(scale(t(yhat_L1_Common)))
  )
  pheatmap(
    combined_mat,
    cluster_cols = FALSE, cluster_rows = TRUE,
    show_rownames = TRUE, show_colnames = FALSE,
    border_color = "white",
    filename = file.path(RESULTS_DIR, "Figure4_track_A.pdf")
  )
}

# ==============================================================================
# 9. Gene Trend Plots (Smoothed Lines)
# ==============================================================================

message("Plotting individual gene trends...")

# Helper function to plot smoothed expression trends
# This reduces code duplication significantly
plot_gene_trend <- function(model, counts_obj, gene, colors = curve_cols) {
  
  p <- plotSmoothers(
    model,
    counts = counts_obj,
    lwd = 0.5,
    curvesCols = colors,
    gene = gene
  ) +
    scale_color_manual(values = colors) +
    ggtitle(gene) +
    theme(
      legend.position = c(0.85, 0.85),
      plot.title = element_text(size = 20, face = "bold", hjust = 0),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16)
    )
  
  return(p)
}

# Define Gene Sets
genes_phase_1 <- c("ANXA1","RAN","MGST3","ACP3","SEC11A","B2M","RPLP1","PPIE")
genes_phase_2 <- c("BCAP31","PRDX3","VAMP2","COPB2","SLC2A14","ARPC3","SND1","DIAPH1","ARAP1","ANXA1","LGALS3","VAMP7","B2M")
genes_phase_3 <- c("DDX39B","PPIE","RPL6","NNT","VAMP8","RPN2","ARHGEF1","RPS27A","SEPTIN6","ADAM10","CNST","DPYSL2","SPTBN1","SCAMP4")

# Figure 4d
# Generate Plots using lapply and the helper function
plots_phase_1 <- lapply(genes_phase_1, function(g) {
  p <- plot_gene_trend(sce_gam, counts(CD4T_sce), g)
  p + ylim(c(0, 0.15)) # Default ylim
})

# Custom Y-limits for specific genes in Phase 1 (as per original code)
# Adjusting specific indices based on gene names
plots_phase_1[[1]] <- plots_phase_1[[1]] + ylim(c(0, 1.25)) # ANXA1
plots_phase_1[[2]] <- plots_phase_1[[2]] + ylim(c(0, 0.25)) # RAN
plots_phase_1[[4]] <- plots_phase_1[[4]] + ylim(c(0, 0.05)) # ACP3
plots_phase_1[[6]] <- plots_phase_1[[6]] + ylim(c(1.9, 2.4)) # B2M
plots_phase_1[[7]] <- plots_phase_1[[7]] + ylim(c(0.9, 1.6)) # RPLP1

# Arrange and Save Phase 1
layout_matrix <- rbind(c(1,2,3,4), c(5,6,7,8))
p_grid_phase1 <- grid.arrange(grobs = plots_phase_1, layout_matrix = layout_matrix)

# ggsave(file.path(RESULTS_DIR, "Fig4d_supp_Phase1.pdf"), p_grid_phase1, width = 18, height = 8)

# Generate Plots for Phase 2 & 3
plots_phase_2 <- lapply(genes_phase_2, function(g) plot_gene_trend(sce_gam, counts(CD4T_sce), g) + ylim(c(0, 0.15)))
plots_phase_3 <- lapply(genes_phase_3, function(g) plot_gene_trend(sce_gam, counts(CD4T_sce), g) + ylim(c(0, 0.15)))

message("Analysis complete. Check the 'results' folder for outputs.")
