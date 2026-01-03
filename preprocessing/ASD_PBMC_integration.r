################################################################################
# Title: Single-cell RNA-seq Preprocessing Pipeline (COVID-19 & Autism)
# Description: 
#   This script integrates reference COVID-19 datasets, in-house Autism/Treat 
#   samples, and public healthy controls (GSE135779). It performs merging, 
#   normalization (SCTransform), and cell cycle scoring.
#
# Requirements: Seurat, reticulate, future
# Author: Jaemyung Jang / Korea Brain Research Instite (piloter2@kbri.re.kr)
################################################################################

# 1. Setup Environment ---------------------------------------------------------
library(Seurat)
library(reticulate)
library(future)

# Configure Python (Adjust path as needed or rely on default environment)
# use_python("/usr/bin/python3") 
py_discover_config()

# Set global options
options(future.globals.maxSize = 32 * 1024 * 1024 * 1024)

# Define Base Paths (Change these paths for your local environment)
DATA_DIR <- "./data" 
PREPROC_DIR <- "./preprocessing"
GSE_DIR <- file.path(PREPROC_DIR, "GSE135779")

# Load custom functions
# Ensure this file is included in your repository or documented
source(file.path(DATA_DIR, "single_cell_function_HUMAN.r")) 


# 2. Utility Functions ---------------------------------------------------------

#' Fix Gene Names
#' Replaces hyphens with dots in row names (Seurat compatibility)
#' @param obj Seurat object
#' @return Seurat object with sanitized gene names
fix_gene_names <- function(obj) {
  # Fix counts dimnames
  counts_names <- rownames(GetAssayData(obj, slot = "counts"))
  new_names <- gsub("\\-", ".", counts_names)
  
  # Update matrix row names (Accessing internals depending on Seurat version)
  # Note: Ideally, fix names on the raw matrix before creating the Seurat Object
  # but adapting to the logic provided:
  if (!is.null(obj@assays$RNA@counts)) {
    rownames(obj@assays$RNA@counts) <- new_names
  }
  if (!is.null(obj@assays$RNA@data)) {
    rownames(obj@assays$RNA@data) <- new_names
  }
  if (!is.null(obj[["RNA"]]@meta.features)) {
    rownames(obj[["RNA"]]@meta.features) <- gsub("\\-", ".", rownames(obj[["RNA"]]@meta.features))
  }
  return(obj)
}

# 3. Load Reference Data (Blish COVID) -----------------------------------------
message("Loading reference dataset...")
reference_dataset <- readRDS(file.path(PREPROC_DIR, "blish_covid.seu.rds"))
reference_dataset <- UpdateSeuratObject(reference_dataset)

# Extract Healthy Adults from Reference
Idents(reference_dataset) <- "Ventilated"
reference_adults_matrix <- GetAssayData(subset(reference_dataset, ident = "Healthy"))
rownames(reference_adults_matrix) <- gsub("\\-", ".", rownames(reference_adults_matrix))

reference_obj <- CreateSeuratObject(reference_adults_matrix, project = "COVID", min.cells = 10)
reference_obj$id_new <- "1. Reference"

# Add metadata from original object (Ensure dimensions match)
# Note: Simple merging of metadata works if cell names are preserved
reference_obj <- AddMetaData(reference_obj, reference_dataset@meta.data)


# 4. Load In-House Data (Autism/Treat) -----------------------------------------
message("Loading in-house datasets...")

# Define file mapping
files_info <- list(
  list(path = "before_human1_HLA.rds", id = "2. Autism", donor = "human1"),
  list(path = "before_human2_HLA.rds", id = "2. Autism", donor = "human2"),
  list(path = "after_human1_HLA.rds",  id = "3. Treat",  donor = "human1"),
  list(path = "after_human2_HLA.rds",  id = "3. Treat",  donor = "human2")
)

# Load and process in a loop
human_pbmc <- lapply(files_info, function(info) {
  obj <- readRDS(file.path(PREPROC_DIR, info$path))
  obj$id_new <- info$id
  obj$Donor <- info$donor
  obj <- fix_gene_names(obj)
  return(obj)
})


# 5. Load Public Data (GSE135779) ----------------------------------------------
message("Loading GSE135779 datasets...")

# Define sample mapping
gse_samples <- list(
  "hc_1"  = list(dir = "JB17010", prefix = "GSM4029907_JB17010_"),
  "hc_2"  = list(dir = "JB17017", prefix = "GSM4029914_JB17017_"),
  "hc_3"  = list(dir = "JB17018", prefix = "GSM4029915_JB17018_"),
  "hc_4"  = list(dir = "JB18069", prefix = "GSM4029922_JB18069_"),
  "hc_5"  = list(dir = "JB18070", prefix = "GSM4029923_JB18070_"),
  "hc_6"  = list(dir = "JB18077", prefix = "GSM4029930_JB18077_"),
  "hc_7"  = list(dir = "JB18078", prefix = "GSM4029931_JB18078_"),
  "hc_8"  = list(dir = "JB18083", prefix = "GSM4029938_JB18083_"),
  "hc_9"  = list(dir = "JB18084", prefix = "GSM4029939_JB18084_"),
  "hc_10" = list(dir = "JB18085", prefix = "GSM4029936_JB18085_"),
  "hc_11" = list(dir = "JB18086", prefix = "GSM4029937_JB18086_")
)

# Load, clean, and create Seurat objects
hc_list <- lapply(names(gse_samples), function(n) {
  sample_info <- gse_samples[[n]]
  path <- file.path(GSE_DIR, sample_info$dir)
  
  # Assuming Read10X_GEO is defined in single_cell_function_HUMAN.r
  raw_data <- Read10X_GEO(paste0(path, "/"), sample.names = sample_info$prefix)
  
  # Access the first element of the list returned by Read10X
  matrix_data <- raw_data[[sample_info$prefix]][[1]]
  
  # Fix rownames
  rownames(matrix_data) <- gsub("\\-", ".", rownames(matrix_data))
  rownames(matrix_data) <- gsub("\\_", ".", rownames(matrix_data))
  
  # Create Object
  obj <- CreateSeuratObject(matrix_data, project = "GSE135779", min.cells = 10)
  obj$id_new <- "0. Reference"
  obj$orig.ident <- "GSE135779"
  obj$Donor <- n
  
  # Run custom processing (assuming processB is in sourced file)
  if (exists("processB")) {
    obj <- processB(obj)
  } else {
    warning("Function 'processB' not found. Skipping custom processing.")
  }
  
  return(obj)
})


# 6. Integration and Normalization ---------------------------------------------
message("Merging datasets...")

# Combine all objects into a single list to merge
# Note: 'PBMC.merges' logic implies merging everything onto reference_obj
combined_list <- c(hc_list, human_pbmc)

PBMC.merged <- merge(
  x = reference_obj,
  y = combined_list,
  project = 'PBMC',
  merge.data = FALSE
)

message("Running SCTransform and CellCycleScoring...")

# Ensure cell cycle genes are available
if (!exists("s.gene")) s.gene <- cc.genes$s.genes
if (!exists("g2m.gene")) g2m.gene <- cc.genes$g2m.genes

# 1st Pass: Normalize to Score Cell Cycle
PBMC.merged <- SCTransform(
  PBMC.merged, 
  assay = 'RNA', 
  new.assay.name = 'SCT', 
  verbose = TRUE
)

PBMC.merged <- CellCycleScoring(
  PBMC.merged,
  s.features = intersect(rownames(PBMC.merged), s.gene),
  g2m.features = intersect(rownames(PBMC.merged), g2m.gene),
  assay = 'SCT',
  set.ident = TRUE
)

# 2nd Pass: Regress out Cell Cycle and Mitochondrial effects
PBMC.merged <- SCTransform(
  PBMC.merged,
  assay = 'RNA',
  new.assay.name = 'SCT',
  vars.to.regress = c('percent.mt', 'nFeature_RNA', 'S.Score', 'G2M.Score'),
  conserve.memory = TRUE,
  verbose = TRUE
)

:gc()
require(harmony)
require(Seurat)

PBMC.mergedd<-PBMC.merged %>%
  RunPCA(verbose = TRUE, features = setdiff(PBMC.features,grep("^MT",PBMC.features, value = TRUE))) %>%  
  RunHarmony(reduction.use = "pca", group.by.vars="Donor", max.iter.harmony = 100, max.iter.cluster=250) %>% #"platform",
  RunUMAP(reduction = "harmony", umap.method = "umap-learn", dims=1:50) %>%
  FindNeighbors(reduction = "harmony", dims=1:50) 

PBMC.mergedd <- PrepSCTFindMarkers(PBMC.mergedd)
PBMC.subset.adults <- subset(PBMC.mergedd, subset = Donor == "H1" | Donor == "H2"  | Donor == "H3"  | Donor == "H4"  | Donor == "H5"  | Donor == "H6" )
PBMC.subset <- subset(PBMC.mergedd, subset = Donor != "H1" & Donor != "H2" & Donor != "H3" & Donor != "H4" & Donor != "H5" & Donor != "H6" )
PBMC.subset <- PrepSCTFindMarkers(PBMC.subset)

saveRDS(PBMC.subset, "./preprocessing/ASD_PBMC_subset_revision.RDS")

message("Preprocessing complete.")
