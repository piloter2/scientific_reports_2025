################################################################################
# Title: Integrated proteomic and single-cell transcriptomic profiling elucidates immunomodulatory effects of l-serine in autism spectrum disorder
# Description: Analysis of PBMC subsets, marker visualization, and gene scoring 
#              (SFARI, TADA, Inflammatory) for ASD vs Healthy Controls.
# Author: Jaemyung Jang / Korea Brain Research institue (piloter2@kbri.re.kr)
# Date: 2026-01-02
# R Version: 4.5.2
################################################################################

# ==============================================================================
# 1. Setup & Libraries
# ==============================================================================
# Data manipulation and visualization
library(Seurat)
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

# Single-cell specific tools
library(Nebulosa) # For density plots
library(UCell)    # For module scoring

# Set seed for reproducibility
set.seed(123)

# ==============================================================================
# 2. Configuration & Data Loading
# ==============================================================================
# NOTE: Please adjust the paths below to match your local environment
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
RDS_PATH <- "/path/to/your/ASD_PBMC_subset_revision.RDS" 

# Load Seurat Object
if (!file.exists(RDS_PATH)) stop("RDS file not found. Check the path.")
PBMC.subset <- readRDS(RDS_PATH)

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
# 4. Figure 2B: UMAP by Group (Healthy vs ASD)
# ==============================================================================
# Define logic for separating groups (Assuming 'id_new' defines the split)
# Note: Adjust logic if the first element of table is not always the Control group
group_ids <- names(table(PBMC.subsets$id_new))
control_id <- group_ids[1]

# Plot 1: Healthy Control
p2b_hc <- DimPlot(subset(PBMC.subsets, subset = id_new == control_id), 
                  group.by = "cell.type.coarse", 
                  label = TRUE, label.size = 4, pt.size = 0.01, raster = FALSE) +
  scale_color_brewer(palette = "Set3") +
  ggtitle("Healthy Control") +
  theme(legend.position = "none", 
        text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold", hjust = 0), 
        axis.line = element_blank(), axis.title = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank()) + NoAxes()

# Plot 2: ASD Patients
p2b_asd <- DimPlot(subset(PBMC.subsets, subset = id_new != control_id), 
                   group.by = "cell.type.coarse", 
                   label = TRUE, label.size = 4, pt.size = 0.01, raster = FALSE) +
  scale_color_brewer(palette = "Set3") +
  ggtitle("ASD Patients") +
  theme(legend.position = "none", 
        text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold", hjust = 0), 
        axis.line = element_blank(), axis.title = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank()) + NoAxes()

# (Optional) Display plots
# p2b_hc + p2b_asd

# ==============================================================================
# 5. Figure 2C: Density Plots of Key Markers (Nebulosa)
# ==============================================================================
# Define markers for cell types
cell_types <- c("CD4 T", "CD8 T", "Mono", "NK", "B", "Platelet")
marker_list <- list(
  c("IL7R", "CCR7"),   # CD4 T
  c("CD8A", "CD8B"),   # CD8 T
  c("CD14", "MS4A7"),  # Mono
  c("GNLY", "NKG7"),   # NK
  c("CD79A", "MS4A1"), # B
  c("PPBP", "PF4")     # Platelet
)
names(marker_list) <- cell_types
flat_markers <- unlist(marker_list)

# Generate density plots loop
p_density_list <- lapply(seq_along(flat_markers), function(x){
  p <- plot_density(as.SingleCellExperiment(PBMC.subsets, assay = "SCT"), 
                    flat_markers[[x]], 
                    reduction = "UMAP", combine = TRUE, pal = "plasma", 
                    adjust = 0.1, size = 0.1) +
    theme(axis.title = element_blank(), axis.line = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = c(0.9, 0.9), 
          legend.background = element_rect(colour = "transparent"))
  
  # Remove density legend title for cleaner look
  update_labels(p, list(colour = " "))
})

# Arrange plots in a grid
layout_matrix <- rbind(c(1,2,3), c(4,5,6), c(7,8,9), c(10,11,12))
# options(repr.plot.width=12, repr.plot.height=16) # For Jupyter environments
p2c <- grid.arrange(grobs = p_density_list, layout_matrix = layout_matrix)

# ==============================================================================
# 6. Figure 2D: Cell Type Composition Analysis
# ==============================================================================
# Calculate proportions
prop_table <- prop.table(table(PBMC.subsets$cell.type.coarse, PBMC.subsets$id_new), margin = 2)
pts <- as.data.frame(prop_table)
pts$Var1 <- as.character(pts$Var1)

# Add label only if frequency > 5%
pts <- pts %>% mutate(label = ifelse(Freq > 0.05, scales::percent(Freq, accuracy = .1), ""))

p2d <- ggplot(pts, aes(x = Var2, y = Freq, fill = reorder(Var1, 1:nrow(pts)), label = label)) +
  geom_col(position = "fill", width = 0.5) +
  geom_text_repel(position = position_fill(vjust = 0.5), 
                  vjust = 0.2, hjust = -0.1, size = 5,
                  max.overlaps = 100, box.padding = unit(0.25, "lines"),
                  direction = "y") + 
  scale_fill_brewer(palette = "Set3") +
  theme_bw(base_size = 15) +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0),
        legend.position = 'bottom', legend.title = element_blank()) +
  coord_flip()

# ==============================================================================
# 7. Figure 2E: Gene Set Scoring (SFARI & TADA)
# ==============================================================================
# Load External Gene Lists
sfari_path <- file.path(DATA_DIR, "SFARI-Gene_genes_10-31-2019release_02-21-2022export.csv")
tada_path <- file.path(DATA_DIR, "TADA102.xlsx")

SFARI <- read.csv(sfari_path)
tada_info <- read_excel(tada_path, col_names = TRUE, na = "NA")

# Filter Gene Sets
gene_tada_ASD <- as.vector(na.omit(tada_info$gene[tada_info$ASD_vs_DDID == "ASD"]))

# SFARI filtering logic: Scores 1, 2 or Syndromic(1), present in data
SFARI_12S <- SFARI[(SFARI$gene.score %in% c(1, 2) | SFARI$syndromic == 1) & 
                   (SFARI$gene.symbol %in% rownames(PBMC.subset)), ]
SFARI_12S <- na.omit(SFARI_12S)

# Compute Module Scores using UCell
obj_dat <- AddModuleScore_UCell(PBMC.subsets, features = list(SFARI_12S$gene.symbol), 
                                name = "SFARIUC_SCT1", assay = "SCT", slot = "data")
obj_dat <- AddModuleScore_UCell(obj_dat, features = list(gene_tada_ASD), 
                                name = "TADAUC_SCT1", assay = "SCT", slot = "data")

# Boxplots for Gene Scores
# Plot 1: SFARI Score
p2e_sfari <- ggboxplot(obj_dat@meta.data, x = "cell.type.coarse", y = "signature_1SFARIUC_SCT1", 
                       color = "id_new", width = 0.5, size = 0.1, 
                       palette = "lancet", outlier.size = 0.005) +
  rotate_x_text(angle = 45) +    
  stat_compare_means(aes(group = id_new), label = "p.signif", label.y = 0.12) + 
  ggtitle("SFARI gene score") +
  theme(axis.title = element_blank())

# Plot 2: TADA Score
p2e_tada <- ggboxplot(obj_dat@meta.data, x = "cell.type.coarse", y = "signature_1TADAUC_SCT1", 
                      color = "id_new", width = 0.5, size = 0.1, 
                      palette = "lancet", outlier.size = 0.005) +
  rotate_x_text(angle = 45) +
  stat_compare_means(aes(group = id_new), label = "p.signif", label.y = 0.15) + 
  ggtitle("TADA gene score") +
  theme(axis.title = element_blank())

# ==============================================================================
# 8. Supplementary Figure 2B: Inflammatory Cytokines & GO Terms
# ==============================================================================
# Define Cytokine Genes
cytokine.genes <- c("IFNG", "IL1B", "IL6", "IL8", "IL10", "IL12A", "IL12B", "IL13", "TNF",
                    "CXCL10", "CCL2", "CCL3", "CCL4", "ICAM1", "SELE", "SELP")

# Load GO Term Genes
go_path <- file.path(DATA_DIR, "GO0006954_inflammatory_response.txt")
inflam.genes <- fread(go_path, header = FALSE) %>% pull(V2) %>% unique()

# Calculate Scores
obj_dat <- AddModuleScore_UCell(obj_dat, features = list(cytokine.genes), 
                                name = "Inflammatory cytokines", assay = "SCT", slot = "data")
obj_dat <- AddModuleScore_UCell(obj_dat, features = list(inflam.genes), 
                                name = "GO0006954", assay = "SCT", slot = "data")

# Define Comparisons for Statistics
# WARNING: Check levels of 'Status' to ensure correct indices (1,3) and (3,2)
status_levels_cd4 <- names(table(obj_dat@meta.data$Status[obj_dat@meta.data$cell.type.fine == "CD4T"]))
my_comparisons <- list(status_levels_cd4[c(1, 3)], status_levels_cd4[c(3, 2)])

p2s <- ggboxplot(obj_dat@meta.data, 
                 x = "cell.type.coarse", y = "signature_1Inflammatory cytokines", 
                 color = "id_new", width = 0.5, size = 0.1, 
                 palette = "lancet", outlier.size = 0.005) +
  stat_compare_means(method = "t.test", 
                     comparisons = my_comparisons,  
                     label = "p.format",
                     paired = FALSE, hide.ns = FALSE,
                     label.y = c(0.4, 0.45)) + 
  ggtitle("ASD - Chronic Inflammatory Response (GO0002544)") +
  rotate_x_text(angle = 45) +    
  theme(axis.title = element_blank())

# End of Script
