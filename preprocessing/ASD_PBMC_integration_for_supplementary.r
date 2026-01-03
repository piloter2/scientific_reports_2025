################################################################################
# Project: Single Cell RNA-seq Integration (COVID-19 Ref + Autism + GSE279452 + GSE206295 + GSE168732)
# Author: Jaemayung Jang / Korea Brain Research Institute (piloter2@kbri.re.kr)
# Date: 2026-01-02
# Description: Integration of reference datasets, in-house samples, and 
#              public healthy controls using Seurat and Harmony.
################################################################################

# ==============================================================================
# 1. Environment Setup & Configuration
# ==============================================================================

# Libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(WGCNA)
  library(hdWGCNA)
  library(scCustomize)
  library(unix)
  library(gprofiler2)
  library(harmony)
  library(reticulate)
})

# System Resources
rlimit_as(4 * 1e12) # Increase memory limit
options(future.globals.maxSize = 32 * 1024 * 1024 * 1024)

# Python Configuration (Adjust path if using a specific env)
# use_python("/usr/bin/python3") 
py_discover_config()

# --- USER CONFIGURATION (EDIT THESE PATHS) ---
BASE_DIR    <- "./data"
REF_DIR     <- file.path(BASE_DIR, "references")
RAW_DIR     <- file.path(BASE_DIR, "raw")
FUNC_PATH   <- "./functions/single_cell_function_HUMAN.r"

# Load custom functions
if(file.exists(FUNC_PATH)) source(FUNC_PATH)

# ==============================================================================
# 2. Helper Functions
# ==============================================================================

#' Convert ENSG IDs to Gene Symbols
#' @param data_object List containing 10x data
#' @return Data object with updated row names
convert_ensg_to_symbol <- function(data_object) {
  for (sample_name in names(data_object)) {
    current_genes <- rownames(data_object[[sample_name]][[1]])
    ensg_indices <- grep("^ENSG", current_genes)
    
    if (length(ensg_indices) > 0) {
      ensg_ids_to_convert <- current_genes[ensg_indices]
      
      # Batch conversion using gprofiler2
      gres <- gconvert(query = ensg_ids_to_convert, 
                       organism = "hsapiens", 
                       target = "ENSG", 
                       filter_na = FALSE)
      
      mapped_symbols <- gres$name
      
      # Handle NAs: keep original ID if no symbol found
      na_indices <- is.na(mapped_symbols) | mapped_symbols == "N/A" | mapped_symbols == ""
      mapped_symbols[na_indices] <- ensg_ids_to_convert[na_indices]
      
      # Ensure uniqueness
      mapped_symbols <- make.unique(mapped_symbols)
      
      # Update rownames
      current_genes[ensg_indices] <- mapped_symbols
      rownames(data_object[[sample_name]][[1]]) <- current_genes
      
      message(paste("Processed:", sample_name, "-", length(ensg_indices), "genes converted."))
    } else {
      message(paste("Processed:", sample_name, "- No ENSG IDs found."))
    }
  }
  return(data_object)
}

#' Standardize Gene Names (Remove dashes)
#' @param seurat_obj Seurat Object
#' @return Seurat Object with sanitized names
sanitize_gene_names <- function(seurat_obj) {
  # Logic to rename rows in count/data slots would go here if not handled at matrix read
  # Currently handling this during matrix loading steps below
  return(seurat_obj)
}

# ==============================================================================
# 3. Load Reference Data (COVID-19)
# ==============================================================================
message("Loading Reference Dataset...")

ref_file <- file.path(REF_DIR, "blish_covid.seu.rds")
if(!file.exists(ref_file)) stop("Reference file not found.")

reference_dataset <- readRDS(ref_file)
reference_dataset <- UpdateSeuratObject(reference_dataset)

# Subset for Healthy Adults
Idents(reference_dataset) <- "Ventilated"
reference.adults <- GetAssayData(subset(reference_dataset, ident = "Healthy"))
rownames(reference.adults) <- gsub("\\-", ".", rownames(reference.adults))

reference.obj <- CreateSeuratObject(reference.adults, project = "COVID", min.cells = 10)
reference.obj$id_new <- "1. Reference"
reference.obj <- AddMetaData(reference.obj, reference_dataset@meta.data)

# ==============================================================================
# 4. Load In-House Data (Autism/Treat)
# ==============================================================================
message("Loading In-House Datasets...")

# Define file metadata to avoid repetitive code
in_house_files <- list(
  list(file="before_human1_HLA.rds", id="2. Autism", donor="human1"),
  list(file="before_human2_HLA.rds", id="2. Autism", donor="human2"),
  list(file="after_human1_HLA.rds",  id="3. Treat",  donor="human1"),
  list(file="after_human2_HLA.rds",  id="3. Treat",  donor="human2")
)

human_pbmc <- lapply(in_house_files, function(info) {
  path <- file.path(RAW_DIR, info$file)
  if(!file.exists(path)) { warning(paste("File missing:", path)); return(NULL) }
  
  obj <- readRDS(path)
  obj$id_new <- info$id
  obj$Donor <- info$donor
  
  # Sanitize gene names (replace - with .)
  # Note: Modifying slots directly is risky in Seurat V5, prefer RenameGenes if available
  # Keeping user logic for V3/V4 compatibility:
  counts <- GetAssayData(obj, slot = "counts")
  rownames(counts) <- gsub("\\-", ".", rownames(counts))
  
  # Re-create object to ensure consistency
  new_obj <- CreateSeuratObject(counts, meta.data = obj@meta.data)
  return(new_obj)
})

# Filter out NULLs if files were missing
human_pbmc <- Filter(Negate(is.null), human_pbmc)

# ==============================================================================
# 5. Load Public Healthy Controls (GEO)
# ==============================================================================
message("Loading GEO Healthy Controls...")

# Define GEO metadata
geo_meta <- data.frame(
  sample_name = c(
    "GSM8571146_PHC01.", "GSM8571147_PHC02.", "GSM8571148_PHC03.", "GSM8571149_PHC05.",
    "GSM6250006_Con1_", "GSM6250007_Con2_", "GSM6250008_Con3_",
    "GSM5160432_H1_health_", "GSM5160434_H2_health_", "GSM5160435_H3_health_"
  ),
  repo = c(rep("GSE279452", 4), rep("GSE206295", 3), rep("GSE168732", 3)),
  path = c(
    rep("GSE279452/", 4), rep("GSE206295/", 3), rep("GSE168732/", 3)
  )
)

hc_list <- list()

for(i in 1:nrow(geo_meta)) {
  path <- file.path(RAW_DIR, geo_meta$path[i])
  s_name <- geo_meta$sample_name[i]
  
  # Read Data
  # Note: Assuming Read10X_GEO is from your custom sourced script
  raw_data <- Read10X_GEO(path, sample.names = s_name)
  
  # Convert IDs
  raw_data <- convert_ensg_to_symbol(raw_data)
  
  # Extract matrix (assuming list structure from Read10X_GEO)
  mat <- raw_data[[s_name]][[1]]
  
  # Clean Rownames
  rownames(mat) <- make.unique(gsub("\\-", ".", rownames(mat)))
  
  # Create Seurat Object
  obj <- CreateSeuratObject(mat, project = geo_meta$repo[i], min.cells = 10)
  obj$id_new <- "0. Reference"
  obj$orig.ident <- geo_meta$repo[i]
  obj$Donor <- paste0(geo_meta$repo[i], "_", i)
  
  hc_list[[i]] <- obj
}

# Run custom processing on HC data
# Assuming processA is defined in single_cell_function_HUMAN.r
if(exists("processA")) {
  hc.processed <- lapply(hc_list, processA)
} else {
  hc.processed <- hc_list
  warning("Function processA not found, using raw objects.")
}

# ==============================================================================
# 6. Merging & Normalization (SCTransform)
# ==============================================================================
message("Merging and Normalizing Data...")

# Merge all datasets
PBMC.merges <- merge(
  x = reference.obj, 
  y = c(hc.processed, human_pbmc), 
  project = 'PBMC', 
  merge.data = FALSE
)

# 1. First Pass SCT (To identify cell cycle phase)
# Note: Ensure s.gene and g2m.gene are loaded from your source file or Seurat
if(!exists("s.gene")) data("cc.genes"); s.gene <- cc.genes$s.genes
if(!exists("g2m.gene")) data("cc.genes"); g2m.gene <- cc.genes$g2m.genes

PBMC.merged <- SCTransform(
  PBMC.merges, 
  assay = 'RNA', 
  new.assay.name = 'SCT', 
  vars.to.regress = c('percent.mt', 'nFeature_RNA'),
  verbose = TRUE
)

# 2. Score Cell Cycle
PBMC.merged <- CellCycleScoring(
  PBMC.merged,
  s.features = intersect(rownames(PBMC.merges), s.gene),
  g2m.features = intersect(rownames(PBMC.merges), g2m.gene),
  assay = 'SCT',
  set.ident = TRUE
)

# 3. Second Pass SCT (Regress out Cell Cycle)
PBMC.merged <- SCTransform(
  PBMC.merged,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'S.Score', 'G2M.Score'),
  conserve.memory = TRUE,
  verbose = TRUE
)

# ==============================================================================
# 7. Integration & Dimensionality Reduction
# ==============================================================================
message("Running Harmony Integration...")

# Feature Selection
integ_features <- SelectIntegrationFeatures(
  object.list = c(list(reference.obj), hc.processed, human_pbmc), 
  nfeatures = 5000
)

# Remove Mito/ENSG genes from features
clean_features <- setdiff(integ_features, grep("^MT|^ENSG", integ_features, value = TRUE))

pcs <- 1:50

# Run Pipeline
PBMC.integrated <- PBMC.merged %>%
  RunPCA(verbose = TRUE, features = clean_features) %>% 
  RunHarmony(
    reduction.use = "pca", 
    group.by = "Donor", 
    assay.use = "SCT", 
    max.iter.harmony = 100,
    max.iter.cluster = 250
  ) %>%
  RunUMAP(
    reduction = "harmony", 
    umap.method = "umap-learn", 
    assay = "SCT", 
    dims = pcs
  ) %>%
  FindNeighbors(reduction = "harmony", assay = "SCT", dims = pcs)

# ==============================================================================
# 8. Clustering & Subsetting
# ==============================================================================
message("Final Clustering...")

# Prepare markers for future analysis
PBMC.integrated <- PrepSCTFindMarkers(PBMC.integrated)

# Subset adults (H1-H6 are presumed exclusions based on original code)
# Logic: Keep everything NOT H1-H6
PBMC.subset <- subset(
  PBMC.integrated, 
  subset = !(Donor %in% c("H1", "H2", "H3", "H4", "H5", "H6"))
)

# Find Clusters
PBMC.subset <- FindClusters(
  PBMC.subset, 
  resolution = c(0.8, 2), 
  algorithm = 2, 
  n.start = 10, 
  n.iter = 100
)

# Save Result
# saveRDS(PBMC.subset, file = "results/PBMC_integrated_clustered.rds")

message("Pipeline Completed Successfully.")
