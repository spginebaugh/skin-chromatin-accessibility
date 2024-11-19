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

library(harmony)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed_data/clustered_annotated_seurat.qs")
signac <- qread("data/processed_data/signac_all_additional_info.qs")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Prepare for integration                         ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(signac) <- "ATAC_RNA"
signac <- NormalizeData(signac)
signac <- ScaleData(signac, features = rownames(signac))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Integrate                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## identify anchors
transfer_anchors <- FindTransferAnchors(
  reference = seurat, 
  query = signac,
  features = VariableFeatures(object = seurat),
  reference.assay = "Corrected",
  query.assay = "ATAC_RNA",
  reduction = "cca",
  normalization.method = "LogNormalize"
  )
  
## annotate scATAC-seq via label transfer
celltype_predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = seurat$annotation_level1,
  weight.reduction = signac[["lsi"]],
  dims = 2:30
)

signac <- AddMetaData(signac, metadata = celltype_predictions)

## impute scATAC-seq based on scRNA-seq
genes_use <- VariableFeatures(seurat)
refdata <- GetAssayData(seurat, assay = "Corrected", slot = "data")[genes_use,]

imputation <- TransferData(
  anchorset = transfer_anchors,
  refdata = refdata,
  weight.reduction = signac[["lsi"]],
  dims = 2:30
  )

signac[["RNA_imputed"]] <- imputation


## coembed for plotting
signac[["Corrected"]] <- imputation
DefaultAssay(signac) <- "Corrected"
coembed <- merge(x = seurat[genes_use, ], y = signac[genes_use,])
coembed <- ScaleData(coembed, features = genes_use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes_use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

coembed$sequence_type <- "RNA"
coembed$sequence_type[colnames(coembed) %in% colnames(signac)] <- "ATAC"

qsave(coembed, "data/processed_data/coembed.qs")


## use coembed to help with labeling 

signac$annotation_level1 <- plyr::revalue(signac$ATAC_snn_res.0.8,
                                          c('0' = "Keratinocyte",
                                            '1' = "Lymphoid",
                                            '2' = "Keratinocyte",
                                            '3' = "Fibroblast",
                                            '4' = "Muscle",
                                            '5' = "Vascular_endo",
                                            '6' = "Fibroblast",
                                            '7' = "Myeloid",
                                            '8' = "Keratinocyte",
                                            '9' = "Keratinocyte",
                                            '10' = "Fibroblast",
                                            '11' = "Lymphoid",
                                            '12' = "Vascular_endo",
                                            '13' = "Lymphoid",
                                            '14' = "Muscle",
                                            '15' = "Myeloid", 
                                            '16' = "Fibroblast",
                                            '17' = "Myeloid",
                                            '18' = "Myeloid",
                                            '19' = "Melanocyte",
                                            '20' = "Keratinocyte",
                                            '21' = "Keratinocyte",
                                            '22' = "Lymphatic_endo",
                                            '23' = "Keratinocyte",
                                            '24' = "Keratinocyte",
                                            '25' = "Plasma_B_cell",
                                            '26' = "Keratinocyte"))

qsave(signac, "data/processed_data/signac_integrated.qs")

