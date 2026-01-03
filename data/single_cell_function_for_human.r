################################################################################
# Title: Single-cell Analysis Utility Functions
# Description: Helper functions for loading data, preprocessing, quality control, 
#              doublet detection, and visualization for scRNA-seq analysis.
# Author: jaemyung / Korea Brain Research Institute
# Date: 2026-01-02
# Dependencies: Seurat, SingleCellExperiment, scDblFinder, tidyverse, etc.
################################################################################

# ==============================================================================
# 1. Environment Setup & Library Loading
# ==============================================================================

# Optimize parallel processing for DelayedArray (used by Monocle3/Bioconductor)
options(DelayedArray.block.size = 1000e6)
DelayedArray:::set_verbose_block_processing(TRUE)

# Suppress warnings for cleaner output
options(warn = -1)
options(stringsAsFactors = FALSE)

# Garbage collection to free up memory
base::gc()

#' Load Required Packages
#' Installs missing packages automatically
load_libraries <- function() {
  # CRAN Packages
  cran_packages <- c(
    "Seurat", "SingleCellExperiment", "dplyr", "tidyverse", "data.table", 
    "scales", "RColorBrewer", "ggplot2", "devtools", "cowplot", "patchwork", 
    "ggpubr", "gridExtra", "ggrepel", "ggbeeswarm", "ggthemes", "grDevices", 
    "extrafont", "openxlsx", "readxl", "magrittr", "pbapply", "reshape2", 
    "stringr", "ggvenn", "stringi", "Matrix", "httr", "utils", "tools"
  )
  
  # Bioconductor Packages
  bioc_packages <- c(
    "EnhancedVolcano", "STRINGdb", "glmGamPoi", "slingshot", "ggstance", 
    "clusterProfiler", "org.Hs.eg.db", "enrichplot", "ComplexHeatmap", 
    "scDblFinder", "multienrichjam", "harmony", "dittoSeq", "limma", "biomaRt"
  )
  
  # Install/Load CRAN
  new_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
  if(length(new_cran)) install.packages(new_cran)
  sapply(cran_packages, suppressPackageStartupMessages(require), character.only = TRUE)
  
  # Install/Load Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
  if(length(new_bioc)) BiocManager::install(new_bioc)
  sapply(bioc_packages, suppressPackageStartupMessages(require), character.only = TRUE)
  
  # GitHub Packages
  if (!requireNamespace("nichenetr", quietly = TRUE)) devtools::install_github("saeyslab/nichenetr")
  require("nichenetr")
}

# Run library loading
load_libraries()

# ==============================================================================
# 2. Global Variables & Reference Paths
# ==============================================================================

# NOTE: Adjust these paths to match your repository structure
BASE_DIR <- getwd()
REF_DIR <- file.path(BASE_DIR, "data")

# Cell Cycle Genes (Regev Lab)
cc_file <- file.path(REF_DIR, "regev_lab_cell_cycle_genes.txt")
if(file.exists(cc_file)){
  genes.cc <- readLines(cc_file)
  s.genes <- genes.cc[1:43]
  g2m.genes <- genes.cc[44:97]
} else {
  warning("Cell cycle gene file not found. Please check REF_DIR.")
}

# Mitochondrial Genes (Manual List for specific needs)
mt.genes.custom <- c("^mt.", "Mtarc2", "Mtch1", "Mterf3", "Mterf4",
                     "Mtfmt", "Mtg1", "Mtg2", "Mtif2",
                     "Mto1", "Mtpap", "Mtrf1l")

# ==============================================================================
# 3. Preprocessing & Quality Control Functions
# ==============================================================================

#' Remove Doublets using scDblFinder
#' @param x Seurat object
#' @return Seurat object with doublets removed
doublet_removal <- function(x){
  require(scDblFinder)
  # Convert to SCE for scDblFinder
  sce <- as.SingleCellExperiment(x)
  sce <- scDblFinder(sce)
  
  # Identify doublets
  doublet_cells <- colnames(sce)[which(sce$scDblFinder.class == "doublet")]
  
  # Subset Seurat object
  if(length(doublet_cells) > 0){
    message(paste("Removing", length(doublet_cells), "doublets"))
    x_filtered <- subset(x, cells = setdiff(colnames(x), doublet_cells))
    return(x_filtered)
  } else {
    message("No doublets found")
    return(x)
  }
}

#' Calculate Mitochondrial and Ribosomal Ratios
#' @param x Seurat object
#' @return Seurat object with percent.mt, percent.rb, percent.hb added to metadata
calculate_qc_metrics <- function(x){
  # Hemoglobin genes (exclude specific ones if needed)
  HB_genes <- grep("^HB", rownames(x), value = TRUE)
  exclude_HB <- c("HBS1L", "HBEGF", "HBP2")
  HB_genes <- setdiff(HB_genes, exclude_HB)
  
  # Mitochondrial genes (Standard Human List)
  MT_symb <- c('MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5',
               'MT-ND6','MT-CYB','MT-CO1','MT-CO2','MT-CO3','MT-ATP6',
               'MT-ATP8','MT-RNR2')
  # Intersect with dataset to avoid warnings
  MT_symb <- intersect(MT_symb, rownames(x))
  
  # Calculate percentages
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-") # Standard regex usually safer
  if(length(MT_symb) > 0 & sum(x[["percent.mt"]]) == 0) {
      # Fallback to specific list if regex fails (e.g., different naming convention)
      x[["percent.mt"]] <- PercentageFeatureSet(x, features = MT_symb) 
  }
  
  x[["percent.rb"]] <- PercentageFeatureSet(x, pattern = "^RP[LS]")
  x[["percent.mrp"]] <- PercentageFeatureSet(x, pattern = "^MRP[LS]")
  x[["percent.hb"]] <- PercentageFeatureSet(x, features = HB_genes)
  
  return(x)
}

#' Standard Preprocessing and Filtering (Method A - Adaptive Thresholds)
#' Filters based on MAD (Median Absolute Deviation)
process_sample_adaptive <- function(run) {
  
  # Calculate metrics
  run <- calculate_qc_metrics(run)
  run$log10GenesPerUMI <- log10(run$nFeature_RNA) / log10(run$nCount_RNA)
  
  # Calculate adaptive thresholds (MAD)
  median_nCount <- median(run$nCount_RNA)
  mad_nCount <- mad(run$nCount_RNA)
  
  median_nFeature <- median(run$nFeature_RNA)
  mad_nFeature <- mad(run$nFeature_RNA)
  
  median_percent_MT <- median(run$percent.mt)
  mad_percent_MT <- mad(run$percent.mt)
  
  thresholds_nCount <- c(0, median_nCount + 5*mad_nCount)
  thresholds_nFeature <- c(0, median_nFeature + 5*mad_nFeature)
  thresholds_percent_MT <- c(0, median_percent_MT + 2*mad_percent_MT)
  
  # Apply Filter
  run <- subset(run, subset = (nFeature_RNA >= 500) & 
                  (percent.mt <= thresholds_percent_MT[2]) & 
                  (log10GenesPerUMI > 0.80) & 
                  (nCount_RNA <= thresholds_nCount[2]) & 
                  (nFeature_RNA <= thresholds_nFeature[2]) )
  
  return(run)
}

#' PBMC Specific Preprocessing
#' Lower nFeature threshold (200) often used for PBMCs
process_pbmc <- function(run) {
  
  run <- calculate_qc_metrics(run)
  run$log10GenesPerUMI <- log10(run$nFeature_RNA) / log10(run$nCount_RNA)
  
  # Adaptive thresholds
  median_nCount <- median(run$nCount_RNA)
  mad_nCount <- mad(run$nCount_RNA)
  median_nFeature <- median(run$nFeature_RNA)
  mad_nFeature <- mad(run$nFeature_RNA)
  median_percent_MT <- median(run$percent.mt)
  mad_percent_MT <- mad(run$percent.mt)
  
  thresholds_nCount <- c(0, median_nCount + 5*mad_nCount)
  thresholds_nFeature <- c(0, median_nFeature + 5*mad_nFeature)
  thresholds_percent_MT <- c(0, median_percent_MT + 2*mad_percent_MT)
  
  # Filter (nFeature >= 200 for PBMC)
  run <- subset(run, subset = (nFeature_RNA >= 200) & 
                  (percent.mt <= thresholds_percent_MT[2]) & 
                  (log10GenesPerUMI > 0.80) & 
                  (nCount_RNA <= thresholds_nCount[2]) & 
                  (nFeature_RNA <= thresholds_nFeature[2]) )
  
  return(run)
}


# ==============================================================================
# 4. Visualization Functions
# ==============================================================================

#' Create Polar Bar Chart for Cell Type Proportions
#' @param object Seurat object
#' @param ident Metadata column for grouping (e.g., 'Donor' or 'Condition')
#' @param expr_clus Clusters to explicitly label
celltype_vol_chart <- function (object, ident, expr_clus) {
  resDF <- list()
  
  # Loop through each group in the identity
  for(classification in names(table(object[[ident]]))) {
    
    # Subset cells safely
    cells.use <- colnames(object)[object[[ident]] == classification]
    sub_obj <- subset(object, cells = cells.use)
    
    df <- as.data.frame(table(Idents(sub_obj)))
    colnames(df) <- c("cluster", "counts")
    df$'class' <- classification
    
    threshold <- .05
    
    # Calculate positions for labels
    df <- df %>% 
      arrange(desc(cluster)) %>%
      mutate(prop = round(counts / sum(df$counts), 4)) %>%
      mutate(lab.pos = cumsum(prop) - 0.5*prop) %>%
      mutate(ypos = ifelse(is.na(lab.pos), prop/2, lab.pos),
             xn = ifelse(prop > threshold, 0, .5))
    
    # Create Plot
    p <- ggplot(data = df, aes(x = 2, y = prop, fill = cluster)) +
      geom_bar(stat = "identity") +
      coord_polar("y", start = 0) +
      geom_text_repel(aes(label = ifelse(cluster %in% expr_clus, 
                                         paste(cluster, ":", scales::percent(prop, accuracy=0.1)), 
                                         paste0(cluster)), 
                          y = lab.pos), 
                      color="black", size=5, nudge_x = df$xn, segment.size = .5) +
      theme_void() +
      xlim(.2, 2.5) +
      ggtitle(paste0(classification, " : ", sum(df$counts), " filtered cells")) + 
      theme(plot.title = element_text(color="black", size=12, face="bold")) +
      NoLegend()
    
    resDF[[classification]] <- p
  }
  return(resDF)
}

#' Spearman Correlation Heatmap of Clusters
#' @param obj.combed Seurat object
celltype_cor_chart <- function (obj.combed) {
  # Calculate average expression per cluster
  av.exp <- AverageExpression(obj.combed)$RNA
  
  # Compute Spearman correlation
  cor.exp <- as.data.frame(cor(av.exp, method="spearman"))
  cor.exp$x <- rownames(cor.exp)
  cor.df <- tidyr::gather(data = cor.exp, y, correlation, -x)
  
  # Plot
  p <- ggplot(cor.df, aes(x = x, y = y, fill = correlation)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "red")) + 
    xlab("Cluster") + ylab("Cluster") +
    ggtitle("Spearman correlation matrix") +
    theme_bw() +
    theme(plot.title = element_text(color="black", size=12, face="bold"),
          axis.title = element_text(color="black", size=10),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(p)
}

#' Calculate Silhouette Scores for Clustering Resolution Optimization
#' @param obj Seurat object
check_silhouette_score <- function(obj){
  require(cluster)
  results <- list()
  
  # Calculate distance matrix once (computationally expensive)
  distance_matrix <- dist(Embeddings(obj[['pca']])[, 1:50])
  
  for(res in seq(0.4, 2.0, 0.1)) {
    cluster_col <- paste0("SCT_snn_res.", res)
    
    if(!cluster_col %in% colnames(obj@meta.data)) {
      message(paste("Resolution", res, "not found in metadata. Skipping."))
      next
    }
    
    clusters <- as.numeric(obj@meta.data[[cluster_col]])
    sil <- silhouette(clusters, dist = distance_matrix)
    
    df <- data.frame(
      'resolution' = res, 
      'mean' = mean(sil[,3]),
      'median' = median(sil[,3]),
      'mad' = mad(sil[,3])
    )
    print(df)
    results[[paste0(res)]] <- df
  }
  return(do.call(rbind, results))
}


# ==============================================================================
# 5. Stress Gene Analysis
# ==============================================================================

#' Identify and Filter Stressed Cells
#' @param obj Seurat object
#' @param stg "go" for Gene Ontology list or "dws" for custom list
#' @param cut Quantile cutoff for stress score (default 0.8)
identify_stress_cells <- function(obj, stg="go", cut=0.8){
  require(tidyverse)
  require(moments)
  
  # Load stress genes
  if (stg == "go"){
    # Update path as needed
    genes_path <- file.path(REF_DIR, "DEG_C2.CGP.M10970.txt")
    if(!file.exists(genes_path)) stop("Stress gene file not found.")
    genes.stress <- read_csv(genes_path)
    genes.stress <- as.list(genes.stress[2:nrow(genes.stress), 1]) # Assuming gene names in col 1
  } else if (stg == "dws"){
    genes.stress <- read_delim("./genesets/genes.deg.Stress.csv", ",", escape_double=FALSE, trim_ws=TRUE)
    genes.stress <- as.list(genes.stress[[1]])
  }
  
  # Subset counts to stress genes
  valid_genes <- intersect(rownames(obj), unlist(genes.stress))
  obj.Stress <- GetAssayData(obj, slot = "counts")[valid_genes, ]
  
  # Filter genes with 0 variance
  obj.Stress <- as.matrix(obj.Stress)
  obj.Stress <- obj.Stress[apply(obj.Stress, 1, var) != 0, ]
  
  # PCA on Stress Genes
  obj.Stress.pca <- prcomp(t(obj.Stress), center=TRUE, scale.=TRUE)
  pca_embeddings <- obj.Stress.pca$x[, 1:2]
  colnames(pca_embeddings) <- paste0("Stress", 1:2)
  
  # Adjust sign based on skewness (ensure stress correlates with positive PC)
  if (skewness(pca_embeddings[,1]) < 0) pca_embeddings[,1] <- -pca_embeddings[,1]
  if (skewness(pca_embeddings[,2]) < 0) pca_embeddings[,2] <- -pca_embeddings[,2]
  
  # Add reduction to Seurat object
  obj[["stress"]] <- CreateDimReducObject(embeddings = pca_embeddings, key = "stress", assay = DefaultAssay(obj))
  
  # Define Cutoff
  cdf <- ecdf(pca_embeddings[,1])
  cut.x <- quantile(pca_embeddings[,1], probs=cut)
  stress_cells <- rownames(pca_embeddings)[pca_embeddings[,1] > cut.x]
  
  # Plotting (Optional validation)
  print(DimPlot(obj, reduction = "stress", pt.size = 0.01) + NoAxes() + ggtitle("Stress Score Distribution"))
  print(VlnPlot(obj, features = c("EGR1", "FOS", "JUN"), group.by = "orig.ident")) # Check known stress markers
  
  # Filter Object
  message(paste("Filtering", length(stress_cells), "stressed cells."))
  obj_filtered <- subset(obj, cells = setdiff(colnames(obj), stress_cells))
  
  return(obj_filtered)
}

# ==============================================================================
# 6. Data Loading Helpers (GEO/10X)
# ==============================================================================

#' Custom 10X Data Reader for GEO Downloads
#' Handles various naming conventions of downloaded GEO files
Read10X_GEO <- function(data.dir = NULL, sample.names = NULL, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE) {
  # ... [Implementation as provided in original code, formatted] ...
  # Note: Standard Seurat::Read10X often handles this now, but kept for legacy support of specific file names.
  
  full.data <- list()
  if (is.null(sample.names)) {
    file.list <- list.files(path = data.dir, pattern = "barcodes.tsv", full.names = FALSE)
    sample.names <- gsub("barcodes.tsv.gz", "", file.list)
  }
  
  for (i in seq_along(sample.names)) {
    prefix <- sample.names[i]
    # Construct paths
    barcode.loc <- file.path(data.dir, paste0(prefix, 'barcodes.tsv.gz'))
    features.loc <- file.path(data.dir, paste0(prefix, 'features.tsv.gz'))
    matrix.loc <- file.path(data.dir, paste0(prefix, 'matrix.mtx.gz'))
    
    # Check for genes.tsv (CellRanger < 3.0)
    if (!file.exists(features.loc)) features.loc <- file.path(data.dir, paste0(prefix, 'genes.tsv.gz'))
    
    if (!file.exists(matrix.loc)) stop("Matrix file missing for: ", prefix)
    
    # Read Data
    data <- Matrix::readMM(matrix.loc)
    cell.names <- readLines(barcode.loc)
    feature.names <- read.delim(features.loc, header = FALSE, stringsAsFactors = FALSE)
    
    colnames(data) <- cell.names
    rownames(data) <- make.unique(feature.names[, gene.column])
    
    full.data[[prefix]] <- data
  }
  
  return(full.data)
}
