# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Load functions and libraries                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(scCDC)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_seurat_file <- "data/processed_data/ambient_removed_seurat.qs"

seurat <- qread("data/processed_data/filtered_seurat.qs")

decontam <- seurat

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            run decontamination                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
decontam <- decontam %>% NormalizeData()
decontam <- decontam %>% FindVariableFeatures()
decontam <- decontam %>% ScaleData()
decontam <- decontam %>% RunPCA()

decontam <- decontam %>% RunUMAP(reduction = "pca", dims = 1:40)

## view umap to decide if batch correction is necessary
DimPlot(decontam, group.by = c("sample", "age_group"))

## batch correct
decontam <- decontam %>% RunHarmony(group.by.vars = "sample", dims.use = 1:40, assay.use = "RNA_ranger")
decontam <- decontam %>% RunUMAP(reduction = "harmony", dims = 1:40)
DimPlot(decontam, group.by = c("sample", "age_group"))

## cluster
decontam <- decontam %>% FindNeighbors(reduction = "harmony", dims = 1:40)
decontam <- decontam %>% FindClusters(resolution = c(0.2, 0.4))

## select resolution
DimPlot(decontam, group.by = c("RNA_ranger_snn_res.0.2", "RNA_ranger_snn_res.0.4"))

## select snn_res.0.4 as resolution
Idents(decontam) <- decontam$RNA_ranger_snn_res.0.4

decontam <- decontam %>% ScaleData(features = rownames(decontam))

gcgs <- ContaminationDetection(decontam)
