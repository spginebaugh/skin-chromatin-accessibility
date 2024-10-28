# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(SingleCellExperiment)
library(tidyverse)
library(harmony)

library(ArchR)
library(mclust)

library(ggrastr)

library(BSgenome.Hsapiens.UCSC.hg38)

source("utils/archr_helpers.R")
source("utils/matrix_helpers.R")
source("utils/misc_helpers.R")
source("utils/plotting_config.R")

set.seed(43648)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rna_proj <- qread("data/processed_data/clustered_annotated_seurat.qs")
atac_proj <- loadArchRProject("data/processed_data/integrated_archR")


# plot1 --------------------------------------------------
featureSets <- list(
  "Keratinocytes" = c("KRT5", "KRT10", "KRT14", "KRT15"), 
  "Fibroblast" = c("THY1", "COL1A1", "COL11A1"), 
  "Lymphoid" = c("CD3D", "CD8A", "CD4","IKZF2", "CCL5"), 
  "Plasma_B_cells" = c("CD79A"), 
  "APCs" = c("CD14", "CD86", "CD163", "CD1A", "CLEC9A", "XCR1"), 
  "Melanocytes" = c("MITF", "SOX10", "MLANA"), 
  "Endothlial" = c("VWF", "PECAM1", "SELE"), 
  "Lymphatic" = c("FLT4", "LYVE1"),  
  "Muscle" = c("TPM1", "TAGLN", "MYL9"), 
  "Pericyte" = c("TRPC6", "CCL19"), 
  "Mast_cells" = c("KIT", "TPSB2", "HPGD"), 
  "HF_surface_markers" = c("ITGB8", "CD200", "SOX9")
)

count_mat <- GetAssayData(object=rna_proj, slot="counts")
avgPctMat <- avgAndPctExpressed(count_mat, rna_proj$annotation_level1, feature_normalize=TRUE, min_pct=0)

subGenes <- featureSets %>% do.call("c",.)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]

avgPctMat$grp <- unlist(rna.NamedClust)[as.character(avgPctMat$grp)]

dotPlot(avgPctMat, xcol = "grp", ycol = "feature", color_col = "avgExpr", size_col = "pctExpr", 
        yorder = unlist(featureSets), xorder = names(featureSets))

plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "annotation_level1", embedding = "UMAP")
