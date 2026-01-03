################################################################################
# Title: Integrated proteomic and single-cell transcriptomic profiling elucidates immunomodulatory effects of l-serine in autism spectrum disorder
# Description: CD4 T
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
library(ggsci)
library(forcats)
library(reshape2)


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
# 4. CD4+ T Cell Subsetting and Visualization
# ==============================================================================

# Subset CD4+ T cells
CD4T <- subset(PBMC.subsets, subset = cell.type.coarse == "CD4+ T")

# Define Colors
cols_ <- RColorBrewer::brewer.pal(n = 10, name = "Paired")

# --- A. UMAP with Highlighting ---
# Define the counts for the title
n_cells <- comma(length(colnames(subset(CD4T, subset = id_new %in% c("3. Treat", "2. Autism")))), format="d")

# Figure 3A
p3a <- DimPlot(PBMC.subsets, 
        cells.highlight = colnames(CD4T),
        sizes.highlight = 0.25,
        cols.highlight = cols_[1],
        reduction = "umap", 
        group.by = "cell.type.coarse", 
        label = TRUE, repel = TRUE, 
        label.size = 4, 
        cols = rep("lightgray", length(unique(PBMC.subsets$cell.type.coarse)))) + 
  theme_void() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0),
        legend.position = "none") + 
  ggtitle(paste0("CD4 T cells : ", n_cells, " cells"))

# Supplementary Figure
ps1 <- DotPlot(CD4T,
            features = c("CCR7","PECAM1","SELL","TUBA1B","CDK6","STMN1",
                         "CXCR3","FAS","IL2RA","GZMA","NKG7","PTPRC","CD27"),
            group.by = "cell.type.fine", 
            cols = c("RdBu"), dot.scale = 8) + 
  RotatedAxis() + coord_flip() +
  labs(y = "Subcluster", x = "Gene symbol")

# Specific Subcluster Markers
ps2 <- DotPlot(subset(CD4T, idents = c("CD4n T-A", "CD4n T-B", "CD4n T-C")), # Assuming idents match
            features = c("IL7R","IL18R1","CD97","CD5","STK38","PIM1","LY6G5C"),
            group.by = "cell.type.fine", 
            cols = c("RdBu"), dot.scale = 8) + 
  RotatedAxis() + coord_flip()

require("scProportionTest")
prop_test <- sc_utils(PBMC.subsets)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "cell.type.fine",
  sample_1 = "1. Reference", 
  sample_2 = "2. Autism",
  sample_identity = "id_new"
)

# Filter results for CD4 only to clean up plot
require("scProportionTest")
prop_test <- sc_utils(PBMC.subset)
prop_test <- permutation_test(
  prop_test, 
  cluster_identity = "cell.type.fine",
  sample_1 = names(table(PBMC.subset$id_new))[1], sample_2 = names(table(PBMC.subset$id_new))[2],
  sample_identity = "id_new"
)
prop_test_mod <- prop_test
prop_test_mod@results$permutation <- prop_test@results$permutation[grep("CD4",prop_test@results$permutation$clusters),]

pne_prop<-permutation_plot(prop_test_mod,log2FD_threshold = 0.7, FDR_threshold = 0.01)

prop_test_new <- sc_utils(PBMC.subset)
prop_test_new <- permutation_test(
  prop_test_new, 
  cluster_identity = "cell.type.fine",
  sample_1 = names(table(PBMC.subset$id_new))[2], sample_2 = names(table(PBMC.subset$id_new))[3],
  sample_identity = "id_new"
)
prop_test_mod2 <- prop_test_new
prop_test_mod2@results$permutation <- prop_test_new@results$permutation[grep("CD4",prop_test_new@results$permutation$clusters),]

pne_prop2<-permutation_plot(prop_test_mod2,log2FD_threshold = 0.7, FDR_threshold = 0.01)

# Figure 3C
p3c <- merge(prop_test_mod2@results$permutation[,c(1:3)], prop_test_mod@results$permutation[,c(1:2)], by = "clusters") %>% 
  melt() %>% 
  mutate(values = round(100*value, 2)) %>% 
  group_by(clusters) %>% 
  # mutate(total = sum(values)) %>% 
  ungroup() %>% 
  ggplot(aes(x=clusters, y=values, fill=variable)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_aaas(palette = "default", alpha = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0)) +
  theme(legend.position = 'bottom') +
  theme(legend.title = element_blank()) + 
  coord_flip()

# Figure 3D
pdf(paste0(BASE_DIR. "/results/Fig3d1_,format(Sys.time(), "%Y-%m-%d %H-%M-%S"),".pdf"), width = 5, height =3)
print(pne_prop)
dev.off()
pdf(paste0(BASE_DIR, "/results/Fig3d2_supp_",format(Sys.time(), "%Y-%m-%d %H-%M-%S"),".pdf"), width = 5, height =3)
print(pne_prop2)
dev.off()

# ==============================================================================
# 5. Differential Expression Analysis
# ==============================================================================

#' Wrapper to run FindMarkers with Error Handling
#' @param obj Seurat object
#' @param id1 Group 1 (Test)
#' @param id2 Group 2 (Control)
run_de_analysis <- function(obj, id1, id2) {
  tryCatch({
    obj <- PrepSCTFindMarkers(obj, assay = "SCT")
    Idents(obj) <- "id_new"
    FindMarkers(obj, assay="SCT", test.use = "MAST", 
                ident.1 = id1, ident.2 = id2, 
                group.by = 'id_new')
  }, error = function(e) {
    message("MAST failed, running without recorrect_umi...")
    FindMarkers(obj, assay="SCT", test.use = "MAST", 
                ident.1 = id1, ident.2 = id2, 
                group.by = 'id_new', recorrect_umi = FALSE)
  })
}

# 1. Treat (AST) vs Autism (ASD)
CD4T.markers.AST_vs_ASD <- lapply(unique(CD4T$cell.type.fine), function(t) {
  message(paste("Processing:", t))
  sub_obj <- subset(CD4T, subset = cell.type.fine == t & id_new != "1. Reference")
  run_de_analysis(sub_obj, id1 = "3. Treat", id2 = "2. Autism")
})
names(CD4T.markers.AST_vs_ASD) <- unique(CD4T$cell.type.fine)

# 2. Autism (ASD) vs Healthy (HC)
CD4T.markers.ASD_vs_HC <- lapply(unique(CD4T$cell.type.fine), function(t) {
  message(paste("Processing:", t))
  sub_obj <- subset(CD4T, subset = cell.type.fine == t & id_new != "3. Treat")
  run_de_analysis(sub_obj, id1 = "2. Autism", id2 = "1. Reference")
})
names(CD4T.markers.ASD_vs_HC) <- unique(CD4T$cell.type.fine)

# ==============================================================================
# 6. DE Summary & Rescue Gene Identification
# ==============================================================================

#' Generate Bar Plot for DEGs
plot_deg_counts <- function(marker_list, title) {
  summary_list <- lapply(names(marker_list), function(cluster_name) {
    df <- marker_list[[cluster_name]]
    sig_genes <- df %>% filter(p_val_adj < 0.01)
    data.frame(Cluster = cluster_name, 
               Up = sum(sig_genes$avg_log2FC > 0), 
               Down = sum(sig_genes$avg_log2FC < 0))
  })
  
  plot_data <- do.call(rbind, summary_list) %>%
    pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count") %>%
    mutate(Count_Plot = ifelse(Direction == "Down", -Count, Count))
  
  max_val <- max(abs(plot_data$Count_Plot)) * 1.15
  
  ggplot(plot_data, aes(x = Cluster, y = Count_Plot, fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = abs(Count), vjust = ifelse(Count_Plot >= 0, -0.5, 1.2)), 
              size = 4, fontface = "bold") +
    scale_y_continuous(labels = abs, limits = c(-max_val, max_val)) + 
    geom_hline(yintercept = 0, color = "black") +
    labs(title = title, subtitle = "Threshold: p_adj < 0.01", y = "Count") +
    theme_minimal() + coord_flip()
}

# Figure 3E
p3e1 <- plot_deg_counts(CD4T.markers.ASD_vs_HC, "ASD vs Healthy")
p3e2 <- plot_deg_counts(CD4T.markers.AST_vs_ASD, "Treat vs ASD")

# --- Identification of Rescued Genes (The Reversal Effect) ---
# Focus on Cluster "CD4n T-A"
cluster_int <- "CD4n T-A"

# Genes UP in ASD, DOWN in Treat (Rescue)
insGenes.1 <- setdiff(
  intersect(
    rownames(subset(CD4T.markers.ASD_vs_HC[[cluster_int]], avg_log2FC > 1 & p_val_adj < 0.01)),
    rownames(subset(CD4T.markers.AST_vs_ASD[[cluster_int]], avg_log2FC < -1 & p_val_adj < 0.01))
  ), 
  c("AC008695.1","AC100801.2") # Manual exclusion list
)

# Genes DOWN in ASD, UP in Treat (Rescue)
insGenes.2 <- intersect(
  rownames(subset(CD4T.markers.ASD_vs_HC[[cluster_int]], avg_log2FC < -1 & p_val_adj < 0.01)),
  rownames(subset(CD4T.markers.AST_vs_ASD[[cluster_int]], avg_log2FC > 1 & p_val_adj < 0.01))
)

# --- Heatmap of Rescued Genes ---
avgexp.CD4T <- AverageExpression(CD4T, group.by = "id_new")$SCT
target_genes <- unique(c(insGenes.1, insGenes.2))
mat <- as.data.frame(avgexp.CD4T[target_genes, ])
mat$gene <- rownames(mat)

# Figure 3F
p3f <- melt(mat) %>% 
  group_by(gene) %>%
  mutate(values = scale(value)) %>% # Z-score scaling by gene
  ggplot(aes(x=variable, y=gene)) + 
  geom_tile(aes(fill = values)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_minimal() +
  labs(title = "Rescue Effect Heatmap", x = "Condition", y = "Gene")

# ==============================================================================
# 7. Pathway Enrichment Analysis (GSEA/ORA)
# ==============================================================================
require(clusterProfiler)
require(ReactomePA)
require(org.Hs.eg.db)

#' Run multiple enrichment analyses on a gene list
#' @param marker_df The output from FindMarkers
run_enrichment_suite <- function(marker_df) {
  
  # Filter and Sort Genes
  selected <- marker_df %>% filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1)
  selected <- selected[!grepl("^MT\\.", rownames(selected)), ]
  
  # Map Symbols to Entrez IDs
  ids <- bitr(rownames(selected), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  # Prepare sorted gene list for GSEA
  gene_list <- selected$avg_log2FC
  names(gene_list) <- rownames(selected)
  gene_list <- gene_list[ids$SYMBOL] # Match order
  names(gene_list) <- ids$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  res_list <- list()
  
  if(length(gene_list) > 0) {
    # 1. GO Terms
    for(ont in c("BP", "MF", "CC")) {
      res_list[[ont]] <- setReadable(
        gseGO(geneList = gene_list, ont = ont, OrgDb = org.Hs.eg.db, 
              keyType = "ENTREZID", pvalueCutoff = 0.5, verbose = FALSE),
        OrgDb = org.Hs.eg.db, keyType = "ENTREZID"
      )
    }
    
    # 2. Reactome
    res_list[["Reactome"]] <- setReadable(
      gsePathway(gene_list, pvalueCutoff = 0.5, verbose = FALSE),
      OrgDb = org.Hs.eg.db, keyType = "ENTREZID"
    )
    
    # 3. WikiPathways
    res_list[["Wiki"]] <- setReadable(
      gseWP(gene_list, organism = "Homo sapiens", pvalueCutoff = 0.5, verbose = FALSE),
      OrgDb = org.Hs.eg.db, keyType = "ENTREZID"
    )
  }
  return(res_list)
}

# Run Analysis for specific cluster
enrich_ASD_HC <- run_enrichment_suite(CD4T.markers.ASD_vs_HC[["CD4n T-A"]])
enrich_Treat_ASD <- run_enrichment_suite(CD4T.markers.AST_vs_ASD[["CD4n T-A"]])

# --- Visualization of Enriched Pathways ---
# Combine results for plotting
combine_enrich_results <- function(res_list) {
  rbind(
    data.frame(res_list[["Reactome"]], group="Reactome"),
    data.frame(res_list[["BP"]], group="GO-BP")
  ) %>% filter(pvalue < 0.05)
}

df_ASD_HC <- combine_enrich_results(enrich_ASD_HC)
df_Treat_ASD <- combine_enrich_results(enrich_Treat_ASD)

# Filter for pathways containing our "Rescued Genes"
target_pattern <- paste0(unique(c(insGenes.1, insGenes.2)), collapse = "|")

plot_pathway_nes <- function(df, title) {
  # Filter pathways that contain the target rescue genes
  df_filt <- df[grep(target_pattern, df$core_enrichment), ]
  
  ggplot(df_filt, aes(NES, fct_reorder(Description, NES), fill = group)) + 
    geom_col() + 
    scale_fill_jama() +
    theme_minimal() + 
    labs(title = title, y = NULL)
}

# Figure 3G
p3g1 <- plot_pathway_nes(df_ASD_HC, "ASD vs Healthy (Pathways with Rescue Genes)")
p3g2 <- plot_pathway_nes(df_Treat_ASD, "Treat vs ASD (Pathways with Rescue Genes)")
