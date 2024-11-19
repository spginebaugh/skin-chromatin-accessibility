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
library(SeuratWrappers)

library(Signac)
library(GenomicRanges)

library(ggrastr)

library(AnnotationHub)
library(BiocParallel)

library(harmony)


library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

library(presto)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Functions                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' creates proper ordering of factors for consistent figures
get_patient_group_factor <- function(patient_group, group_order) {
  patient_group <- factor(patient_group, levels = group_order, ordered = FALSE)
  return(patient_group)
}

#' creates a vector of colors using meta_data dataframe
#' enables consistent coloring even is different DFs have a subset of donors
get_donor_colors <- function(meta_data) {
  axis_colors <- sapply(levels(meta_data$sample_ID), function(x) {
    if (grepl("^C_SD", x)) {
      return(ctr_sd_col)
    } else if (grepl("^C_PB", x)) {
      return(ctr_pb_col)
    } else if (grepl("^AA", x)) {
      return(aa_col)
    }
  })
  return(axis_colors)
}

#' Preps seurat objects for this particular dataset
#' sets the correct ordering for patient_group and sample_ID so figures are consistent
prep_data <- function(seurat_obj) {
  sample_ids <- c(
    paste0("C_SD", 1:7),
    paste0("C_PB", 1:3),
    paste0("AA", 1:8)
  )
  sample_ids <- sample_ids[sample_ids %in% seurat_obj$sample_ID]
  
  seurat_obj$sample_ID <- factor(seurat_obj$sample_ID,
                                 levels = sample_ids,
                                 ordered = FALSE
  )
  return(seurat_obj)
}

cluster_seurat <- function(seurat_obj) {
  seurat_obj <- seurat_obj %>% NormalizeData()
  seurat_obj <- seurat_obj %>% FindVariableFeatures()
  seurat_obj <- seurat_obj %>% ScaleData()
  seurat_obj <- seurat_obj %>% RunPCA()
  seurat_obj <- seurat_obj %>% RunHarmony(group.by.vars = "sample_ID", dims.use = 1:30)
  seurat_obj <- seurat_obj %>% RunUMAP(reduction = "harmony", dims = 1:30)
  
  seurat_obj <- seurat_obj %>% FindNeighbors(reduction = "harmony", dims = 1:30)
  seurat_obj <- seurat_obj %>% FindClusters(resolution = c(0.2,0.4))
  
  return(seurat_obj)
}

cluster_signac <- function(signac_obj){
  DefaultAssay(kera_atac) <- "ATAC"
  
  signac_obj <- RunTFIDF(signac_obj)
  signac_obj <- FindTopFeatures(signac_obj, min.cutoff = 'q0')
  signac_obj <- RunSVD(signac_obj)
  signac_obj <- RunHarmony(signac_obj, group.by.vars = "sample_ID", reduction.use = 'lsi', dims.use = 2:30, assay.use = "ATAC", project.dim = FALSE)
  signac_obj <- RunUMAP(object = signac_obj, reduction = 'harmony', dims = 1:29)
  signac_obj <- FindNeighbors(object = signac_obj, reduction = 'harmony', dims = 1:29)
  signac_obj <- FindClusters(object = signac_obj, 
                            algorithm = 3, 
                            resolution = c(0.6,0.8,1))
  return(signac_obj)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed_data/clustered_annotated_seurat.qs")
signac <- qread("data/processed_data/signac_integrated.qs")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                organize data                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat$patient_group <- plyr::revalue(seurat$disease_status, 
                                      c("control surgical tissue" = "C_SD",
                                        "control punch biopsy" = "C_PB",
                                        "alopecia areata punch biopsy" = "AA"))

signac$patient_group <- plyr::revalue(signac$disease_status, 
                                      c("control surgical tissue" = "C_SD",
                                        "control punch biopsy" = "C_PB",
                                        "alopecia areata punch biopsy" = "AA"))


# fix signac sample columns so they are consistent with seurat
tmp_sample <- signac$sample
tmp_sampleid <- signac$sample_ID
signac$sample_ID <- tmp_sample
signac$sample <- tmp_sampleid


group_order <- c("C_SD","C_PB","AA")
seurat$patient_group <- get_patient_group_factor(seurat$patient_group, group_order)
signac$patient_group <- get_patient_group_factor(signac$patient_group, group_order)

signac@meta.data <- signac@meta.data %>% dplyr::rename(
  ms_broad_ct = manuscript_BroadClust,
  ms_fine_ct = manuscript_NamedClust,
  ms_veryfine_ct = manuscript_FineClust
)

# set colors for plotting
ctr_sd_col <- "#1a53ff"
ctr_pb_col <- "#7ea715"
aa_col <- "#7c1158"

colors <- c(
  "C_SD" = ctr_sd_col,
  "C_PB" = ctr_pb_col,
  "AA" = aa_col
)

seurat <- prep_data(seurat)
signac <- prep_data(signac)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                subclustering                             ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## keratinocytes -------------------------
  ### RNA
kera_rna <- seurat[,seurat$annotation_level1 %in% c("Keratinocyte", "Trichocyte")] 

kera_rna <- cluster_seurat(kera_rna)

DimPlot(kera_rna, group.by = c("sample_ID", "patient_group"))
DimPlot(kera_rna, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(kera_rna, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"), raster = FALSE, label = TRUE)
FeaturePlot(kera_rna, min.cutoff = "q1", max.cutoff = "q99", 
            features = c("cxds_score","bcds_score","hybrid_score","DF_score"))

DimPlot(kera_rna, group.by = c("Corrected_snn_res.0.4"), raster = FALSE, label = FALSE,
        split.by = "patient_group")

kera_rna <- kera_rna[,!(kera_rna$Corrected_snn_res.0.2 %in% c(9,10))]
kera_rna <- cluster_seurat(kera_rna)

  ### ATAC
kera_atac <- signac[,signac$annotation_level1 %in% c("Keratinocyte", "Trichocyte")] 
DefaultAssay(kera_atac) <- "ATAC"

kera_atac <- cluster_signac(kera_atac)
DimPlot(kera_atac, group.by = c("sample_ID", "patient_group"))
DimPlot(kera_atac, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(kera_atac, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"), raster = FALSE, label = TRUE)




