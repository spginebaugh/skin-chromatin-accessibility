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

library(Signac)
library(GenomicRanges)

library(ggrastr)

library(AnnotationHub)
library(BiocParallel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed_data/signac_merge.qs")

dblt <- qread("data/processed_data/amulet_doublets.qs")

meta_sample <- read_csv("manuscript_metadata/manuscript_table_S1_sample_info.csv")
meta_cell_atac <- read_csv("manuscript_metadata/manuscript_table_S2_scatac_meta.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Organize                                ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in names(dblt)){
  dblt[[i]]$sample <- i
  dblt[[i]]$barcode <- paste0(i, "_", rownames(dblt[[i]]))
}
dblt <- do.call(rbind, dblt)
rownames(dblt) <- dblt$barcode
colnames(dblt)[1:6] <- paste0("amulet_", colnames(dblt)[1:6])
