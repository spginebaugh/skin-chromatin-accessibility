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
decontam <- decontam %>% FindClusters(resolution = c(0.2))
DimPlot(decontam, group.by = c("RNA_snn_res.0.2"), label = TRUE)
## select snn_res.0.2 as resolution for clusters
Idents(decontam) <- decontam$RNA_snn_res.0.2



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            run decontamination                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## convert seurat to older version to work with scCDC
options(Seurat.object.assay.version = "v3")
decontam_v3 <- CreateSeuratObject(counts = GetAssayData(decontam, slot = "counts", assay = "RNA"),
                                  meta.data = decontam@meta.data)

## rename clusters to work with program
decontam_v3$seurat_clusters <- factor(paste0("g",decontam_v3$RNA_snn_res.0.2))
decontam_v3@active.ident <- ((decontam_v3$seurat_clusters))

## run decontamination
gcgs <- ContaminationDetection(decontam_v3)
cont_ratio <- ContaminationQuantification(decontam_v3, rownames(gcgs))
seuratobj_corrected = ContaminationCorrection(decontam_v3,rownames(gcgs))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     add info to original seurat object                   ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options(Seurat.object.assay.version = "v5")
seurat <- seurat[, colnames(seurat) %in% colnames(seuratobj_corrected)]

seurat[["Corrected"]] <- CreateAssay5Object(counts = GetAssayData(seuratobj_corrected, slot ="counts", assay = "Corrected"))
DefaultAssay(seurat) <- "Corrected"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    save                                  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave(seurat, output_seurat_file)
