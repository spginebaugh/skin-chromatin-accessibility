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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_seurat_file <- "data/processed_data/clustered_annotated_seurat.qs"
load("data/processed_data/cycle.rda")

seurat <- qread("data/processed_data/ambient_removed_seurat.qs")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  cluster                                 ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- seurat %>% NormalizeData()
seurat <- seurat %>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes)
seurat <- seurat %>% FindVariableFeatures()
seurat <- seurat %>% ScaleData()
seurat <- seurat %>% RunPCA()
seurat <- seurat %>% RunHarmony(group.by.vars = "sample", dims.use = 1:40)
seurat <- seurat %>% RunUMAP(reduction = "harmony", dims = 1:40)

DimPlot(seurat, group.by = c("sample"))
DimPlot(seurat, group.by = c("Phase", "sex","disease_status"))
DimPlot(seurat, group.by = c("ms_broad_ct", "ms_fine_ct"))

seurat <- seurat %>% FindNeighbors(reduction = "harmony", dims = 1:40)
seurat <- seurat %>% FindClusters(resolution = c(0.2,0.4))

DimPlot(seurat, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"))

## selected resolution of 0.2
Idents(seurat) <- seurat$Corrected_snn_res.0.2
DimPlot(seurat, group.by = c("Corrected_snn_res.0.2","ms_broad_ct"), label = TRUE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  annotate                                ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat$annotation_level1 <- plyr::revalue(seurat$Corrected_snn_res.0.2, 
                                    c('0' = "",
                                      '1' = "",
                                      '2' = "",
                                      '3' = "",
                                      '4' = "",
                                      '5' = "",
                                      '6' = "",
                                      '7' = "",
                                      '8' = "",
                                      '9' = "",
                                      '10' = "",
                                      '11' = "",
                                      '12' = "",
                                      '13' = "",
                                      '14' = "",
                                      '15' = ""))
