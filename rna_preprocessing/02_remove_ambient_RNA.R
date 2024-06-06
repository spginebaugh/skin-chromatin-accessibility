# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Load functions and libraries                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(scCDC)
library(harmony)
set.seed(43648)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_seurat_file <- "data/processed_data/ambient_removed_seurat.qs"

seurat <- qread("data/processed_data/filtered_seurat.qs")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                preprocessing                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter out low expression genes
seurat <- seurat[rowSums(GetAssayData(seurat, slot = "counts", assay = "RNA") > 1) > 10,]

decontam <- seurat

decontam <- decontam %>% NormalizeData()
decontam <- decontam %>% FindVariableFeatures()
decontam <- decontam %>% ScaleData()
decontam <- decontam %>% RunPCA()

## batch correct
decontam <- decontam %>% RunHarmony(group.by.vars = "sample", dims.use = 1:40)
decontam <- decontam %>% RunUMAP(reduction = "harmony", dims = 1:40)
DimPlot(decontam, group.by = c("sample", "sex","disease_status"))
DimPlot(decontam, group.by = c("ms_broad_ct", "ms_fine_ct"))
DimPlot(decontam, group.by = c("scDbl_class", "DF_classification"))

## cluster
decontam <- decontam %>% FindNeighbors(reduction = "harmony", dims = 1:40)
decontam <- decontam %>% FindClusters(resolution = c(0.6))

## select resolution
DimPlot(decontam, group.by = c("RNA_snn_res.0.6"), label = TRUE)

## remove doublets 
## based on separation of doublet clusters, 
## we are going to be aggressive in doublet removal
decontam <- decontam[,!(decontam$RNA_snn_res.0.6 %in% c(14,24))]
decontam <- decontam[, (decontam$scDbl_class == "singlet") & (decontam$DF_classification == "Singlet")]

## recluster for input to scCDC
decontam <- decontam %>% FindVariableFeatures()
decontam <- decontam %>% ScaleData()
decontam <- decontam %>% RunPCA()

## batch correct
decontam <- decontam %>% RunHarmony(group.by.vars = "sample", dims.use = 1:40)
decontam <- decontam %>% RunUMAP(reduction = "harmony", dims = 1:40)
DimPlot(decontam, group.by = c("sample", "sex","disease_status"))
DimPlot(decontam, group.by = c("ms_broad_ct", "ms_fine_ct"))

## cluster
decontam <- decontam %>% FindNeighbors(reduction = "harmony", dims = 1:40)
decontam <- decontam %>% FindClusters(resolution = c(0.2,0.4))
DimPlot(decontam, group.by = c("RNA_snn_res.0.2","RNA_snn_res.0.4"), label = TRUE)
## select snn_res.0.2 as resolution for clusters
Idents(decontam) <- decontam$RNA_snn_res.0.2

## scale all genes, as recommended by the scCDC vignette
decontam <- decontam %>% ScaleData(features = rownames(decontam))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            run decontamination                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gcgs <- ContaminationDetection(decontam)
