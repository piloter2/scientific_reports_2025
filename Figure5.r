################################################################################
# Title: CD4+ T Cell Network Analysis (hdWGCNA) & Global DE Profiling
# Description: 
#   1. Export count matrices and metadata for external validation.
#   2. Perform Differential Expression (DE) analysis (AST vs ASD, ASD vs HC).
#   3. Construct Co-expression Networks for CD8+ T cells using hdWGCNA.
# Author: Jaemyung Jang / Korea Brain Research Institute (piloter2@kbri.re.kr)
# Date: 2026-01-02
# R Version: 4.5.2
################################################################################

# ==============================================================================
# 1. Setup & Configuration
# ==============================================================================

# --- Libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  
  # WGCNA related
  library(WGCNA)
  library(hdWGCNA)
  
  # Enrichment
  library(enrichR)
  library(GeneOverlap)
})

# --- Parameters & Paths ---
BASE_DIR    <- getwd() 
DATA_DIR <- file.path(BASE_DIR, "data")
RDS_PATH <- "/path/to/your/ASD_PBMC_subset_revision.RDS" 
OUTPUT_DIR  <- file.path(BASE_DIR, "results/v5.0.0")

# Load Seurat Object
if (!file.exists(RDS_PATH)) stop("RDS file not found. Check the path.")
PBMC.subset <- readRDS(RDS_PATH)

# Ensure directories exist
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# --- Group Definitions (Critical for reproducibility) ---
# NOTE: Adjust these names based on your actual metadata levels
# Using variables avoids errors from changing factor orders
GROUP_CTRL  <- "1. Control"  # Replaces names(table(...))[1]
GROUP_ASD   <- "2. Autism"   # Replaces names(table(...))[2]
GROUP_TREAT <- "3. Treat"    # Replaces names(table(...))[3]

# ==============================================================================
# 2. Export Data for External Use
# ==============================================================================
message("Exporting count matrices and metadata...")

# Helper to export counts
export_counts <- function(seurat_obj, group_id, filename) {
  subset_obj <- subset(seurat_obj, subset = id_new == group_id)
  counts <- as.data.frame(GetAssayData(subset_obj, layer = "counts", assay = "SCT"))
  write.table(counts, file.path(OUTPUT_DIR, filename), sep='\t', quote=FALSE)
}

# Helper to export metadata
export_meta <- function(seurat_obj, group_id, filename) {
  subset_obj <- subset(seurat_obj, subset = id_new == group_id)
  meta_df <- data.frame(
    "barcode_sample" = colnames(subset_obj),
    "cell_type"      = subset_obj@meta.data$cell.type.coarse
  )
  write.table(meta_df, file.path(OUTPUT_DIR, filename), 
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Execute Exports (Treat/AO and ASD/BO)
# AO: After Om (Treatment), BO: Before Om (Autism)
export_counts(PBMC.subset, GROUP_TREAT, "ASD_PBMC_final_supp_PBMC_AO_count.txt")
export_counts(PBMC.subset, GROUP_ASD,   "ASD_PBMC_final_supp_PBMC_BO_count.txt")
export_meta(PBMC.subset, GROUP_ASD,   "ASD_PBMC_final_supp_PBMC_BO_meta.tsv")
export_meta(PBMC.subset, GROUP_TREAT, "ASD_PBMC_final_supp_PBMC_AO_meta.tsv")

# ==============================================================================
# 3. Differential Expression Analysis (Robust Wrapper)
# ==============================================================================
message("Running Differential Expression Analysis...")

# Robust DE function with retry logic for MAST errors
perform_de_with_retry <- function(seurat_obj, cell_type, ident1, ident2, group_col = "id_new") {
  
  # Subset for specific cell type
  obj_sub <- subset(seurat_obj, subset = cell.type.coarse == cell_type)
  
  # Check cell counts
  n_g1 <- sum(obj_sub[[group_col]] == ident1)
  n_g2 <- sum(obj_sub[[group_col]] == ident2)
  
  if (n_g1 < 3 | n_g2 < 3) {
    message(sprintf("Skip %s: Not enough cells (%s: %d, %s: %d)", cell_type, ident1, n_g1, ident2, n_g2))
    return(NULL)
  }
  
  # Further subset to remove irrelevant groups for cleaner processing
  obj_sub <- subset(obj_sub, subset = !!sym(group_col) %in% c(ident1, ident2))
  
  tryCatch({
    obj_sub <- PrepSCTFindMarkers(obj_sub, assay = "SCT")
    Idents(obj_sub) <- group_col
    
    # Attempt 1: Standard MAST
    FindMarkers(obj_sub, assay = "SCT", test.use = "MAST",
                ident.1 = ident1, ident.2 = ident2, group.by = group_col)
    
  }, error = function(e) {
    message(sprintf("Error in %s: %s", cell_type, e$message))
    tryCatch({
      message(" -> Retrying with recorrect_umi = FALSE...")
      # Attempt 2: MAST with recorrect_umi = FALSE
      FindMarkers(obj_sub, assay = "SCT", test.use = "MAST",
                  ident.1 = ident1, ident.2 = ident2, group.by = group_col,
                  recorrect_umi = FALSE)
    }, error = function(e2) {
      message(" -> Retry failed.")
      return(NULL)
    })
  })
}

# Run DE: Treatment (AST) vs Autism (ASD)
cell_types <- names(table(PBMC.subset$cell.type.coarse))
PBMC.markers.AST_vs_ASD <- lapply(cell_types, function(ct) {
  perform_de_with_retry(PBMC.subset, ct, GROUP_TREAT, GROUP_ASD)
})
names(PBMC.markers.AST_vs_ASD) <- cell_types

# Run DE: Autism (ASD) vs Healthy Control (HC)
PBMC.markers.ASD_vs_HC <- lapply(cell_types, function(ct) {
  perform_de_with_retry(PBMC.subset, ct, GROUP_ASD, GROUP_CTRL)
})
names(PBMC.markers.ASD_vs_HC) <- cell_types

# ==============================================================================
# 4. hdWGCNA Analysis for CD8+ T Cells
# ==============================================================================
message("Starting hdWGCNA for CD8+ T cells...")

# --- 4.1 Preprocessing ---
cd8_obj <- subset(PBMC.subset, subset = cell.type.coarse == "CD8+ T" & id_new != "0. Reference")

cd8_obj <- SetupForWGCNA(
  cd8_obj,
  features = VariableFeatures(cd8_obj),
  wgcna_name = "SCT"
)

# Construct Metacells (Aggregation)
cd8_obj <- MetacellsByGroups(
  seurat_obj = cd8_obj,
  group.by = "id_new",
  k = 25,
  max_shared = 12,
  min_cells = 50,
  reduction = 'umap',
  ident.group = 'id_new',
  slot = 'scale.data',
  assay = 'SCT'
)

cd8_obj <- SetDatExpr(cd8_obj)

# --- 4.2 Network Construction ---
# Soft Power Thresholding
cd8_obj <- TestSoftPowers(cd8_obj)
# plot_list <- PlotSoftPowers(cd8_obj)
# ggsave(file.path(OUTPUT_DIR, "softpower.png"), wrap_plots(plot_list, ncol=2))

# Construct Network
cd8_obj <- ConstructNetwork(
  cd8_obj,
  soft_power = 7,  # Selected based on scale-free topology
  tom_name = "SCT_cells",
  overwrite_tom = TRUE
)

# Module Identification
cd8_obj <- ModuleEigengenes(cd8_obj)
cd8_obj <- ModuleConnectivity(cd8_obj)

# Dendrogram Plot
pdf(file.path(CACHE_DIR, "Figure5_CD8_dendrogram.pdf"), width = 10, height = 10)
PlotDendrogram(cd8_obj, main='hdWGCNA SCT Dendrogram')
dev.off()

# Extract Module Info
modules <- GetModules(cd8_obj) %>% subset(module != 'grey')
MEs <- GetMEs(cd8_obj, harmonized=TRUE)
cd8_obj@meta.data <- cbind(cd8_obj@meta.data, MEs)

# --- 4.3 Differential Module Expression (DME) ---
# Comparison: ASD vs Control (within CD8 T cells)
group1_barcodes <- rownames(subset(cd8_obj@meta.data, id_new == GROUP_CTRL))
group2_barcodes <- rownames(subset(cd8_obj@meta.data, id_new == GROUP_ASD))

DMEs <- FindDMEs(
  cd8_obj,
  barcodes1 = group2_barcodes, # ASD
  barcodes2 = group1_barcodes, # Control
  test.use = 'wilcox'
)

# DME Volcano Plot
p_volcano <- DMEs %>%
  ggplot(aes(x=avg_log2FC, y= -log10(p_val_adj), color = module)) +
  geom_point(size = 5) + # Adjusted size for clarity
  xlab("Average Log2 Fold Change (ASD vs Ctrl)") +
  ylab("-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "red", linetype ="dashed") + # Adjusted thresholds
  geom_hline(yintercept = -log10(0.05), colour = "red", linetype ="dashed") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 14)
  )

# --- 4.4 Module Visualization (UMAP) ---
cd8_obj <- RunModuleUMAP(
  cd8_obj,
  n_hubs = 16,
  n_neighbors = 15,
  min_dist = 0.1
)

umap_df <- GetModuleUMAP(cd8_obj)

# Define genes of interest for labeling
ev_genes <- c() # TODO: Load EVgene.pos/neg if available
genes_of_interest <- c(ev_genes, "CD74", "CD8A", "CD8B", "PAG1")
sig_modules <- DMEs %>% dplyr::filter(p_val_adj < 0.05) %>% pull(module)

p_umap <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(color=color, size=kME), alpha=0.8) +
  scale_size_continuous(range = c(0.5, 3)) + 
  scale_color_identity() + 
  # Label specific genes in significant modules
  geom_text_repel(
    data = subset(umap_df, gene %in% genes_of_interest & module %in% sig_modules),
    aes(label = gene),
    max.overlaps = 20
  ) +
  theme_minimal() +
  ggtitle("Module UMAP")

# Hub Gene Network
HubGeneNetworkPlot(
  cd8_obj,
  n_hubs = 16, 
  n_other = 20,
  edge_prop = 0.75,
  mods = sig_modules
)

# ==============================================================================
# 5. Enrichment Analysis (EnrichR)
# ==============================================================================
message("Running Enrichment Analysis...")

dbs <- c('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023')
cd8_obj <- RunEnrichr(cd8_obj, dbs = dbs, max_genes = Inf)

enrich_df <- GetEnrichrTable(cd8_obj)

# Extract gene count from "Overlap" column (e.g., "5/200" -> 5)
enrich_df$counts <- as.numeric(str_split(enrich_df$Overlap, "/", simplify = TRUE)[,1])

# Filter and Plot
target_genes <- c("CD8A", "CD8B", "PAG1", "MTPN", "SPARC", "VAMP7", "SEPTIN6") # Custom list

# Filter for terms containing target genes and belonging to significant modules
filtered_enrich <- enrich_df %>%
  filter(
    grepl(paste(target_genes, collapse = "|"), Genes),
    module %in% sig_modules,
    P.value < 0.01,
    counts > 3,
    db == "GO_Biological_Process_2023"
  )

if(nrow(filtered_enrich) > 0){
  p_enrich <- ggplot(filtered_enrich, aes(x = module, y = forcats::fct_reorder(Term, module))) + 
    geom_point(aes(color = Adjusted.P.value, size = counts)) +
    scale_color_viridis_c(option = "inferno", direction = -1) +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + 
    labs(x = "Co-expression Modules", y = NULL, title = "GO Enrichment (Biological Process)") +
    theme(axis.text.y = element_text(size = 10))
  
} else {
  message("No enrichment terms met the criteria for plotting.")
}

message("Analysis Complete.")
