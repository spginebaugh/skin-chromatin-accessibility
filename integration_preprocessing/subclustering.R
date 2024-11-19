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

# library(cicero)
# library(monocle3)
set.seed(1487)
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
  set.seed(1487)
  seurat_obj <- seurat_obj %>% RunHarmony(group.by.vars = "sample_ID", dims.use = 1:20)
  seurat_obj <- seurat_obj %>% RunUMAP(reduction = "harmony", dims = 1:20)
  set.seed(1487)
  seurat_obj <- seurat_obj %>% FindNeighbors(reduction = "harmony", dims = 1:20)
  seurat_obj <- seurat_obj %>% FindClusters(resolution = c(0.2,0.4))
  
  return(seurat_obj)
}

cluster_signac <- function(signac_obj){
  DefaultAssay(signac_obj) <- "ATAC"
  
  signac_obj <- RunTFIDF(signac_obj)
  signac_obj <- FindTopFeatures(signac_obj, min.cutoff = 'q0')
  signac_obj <- RunSVD(signac_obj)
  set.seed(1487)
  signac_obj <- RunHarmony(signac_obj, group.by.vars = "sample_ID", reduction.use = 'lsi', dims.use = 2:13, assay.use = "ATAC", project.dim = FALSE)
  signac_obj <- RunUMAP(object = signac_obj, reduction = 'harmony', dims = 1:12)
  set.seed(1487)
  signac_obj <- FindNeighbors(object = signac_obj, reduction = 'harmony', dims = 1:12)
  signac_obj <- FindClusters(object = signac_obj, 
                            algorithm = 3, 
                            resolution = c(0.6,0.8,1))
  return(signac_obj)
}

# # decided to exclude from analysis because it is extremely slow
# find_coaccessable_networks <- function(signac_obj){
#   DefaultAssay(signac_obj) <- "ATAC"
#   signac_cds <- as.cell_data_set(signac_obj)
#   signac_cicero <- make_cicero_cds(signac_cds, reduced_coordinates = reducedDims(signac_cds)$UMAP)
#   
#   genome <- seqlengths(signac_obj)
#   genome.df <- data.frame("chr" = names(genome), "length" = genome)
#   conns <- run_cicero(signac_cicero, genomic_coords = genome.df, sample_num = 100)
#   ccans <- generate_ccans(conns)
#   links <- ConnectionsToLinks(conns = conns, ccans = ccans)
#   Links(signac_obj) <- links
#   return(signac_obj)
# }

integrate_and_label_transfer <- function(seurat_obj, signac_obj){
  DefaultAssay(signac_obj) <- "ATAC_RNA"
  signac_obj <- NormalizeData(signac_obj)
  signac_obj <- ScaleData(signac_obj, features = rownames(signac_obj))
  
  transfer_anchors <- FindTransferAnchors(
    reference = seurat_obj, 
    query = signac_obj,
    features = VariableFeatures(object = seurat_obj),
    reference.assay = "Corrected",
    query.assay = "ATAC_RNA",
    reduction = "cca",
    normalization.method = "LogNormalize"
  )
  
  celltype_predictions_clusters <- TransferData(
    anchorset = transfer_anchors,
    refdata = seurat_obj$Corrected_snn_res.0.2,
    weight.reduction = signac_obj[["lsi"]],
    dims = 2:13
  )
  colnames(celltype_predictions_clusters)[1] <- "cluster_transfer"
  signac_obj <- AddMetaData(signac_obj, metadata = celltype_predictions_clusters)
  
  celltype_predictions_ms <- TransferData(
    anchorset = transfer_anchors,
    refdata = seurat_obj$ms_veryfine_ct,
    weight.reduction = signac_obj[["lsi"]],
    dims = 2:13
  )
  colnames(celltype_predictions_ms)[1] <- "ms_veryfine_transfer"
  signac_obj <- AddMetaData(signac_obj, metadata = celltype_predictions_ms)
  
  ## impute scATAC-seq based on scRNA-seq
  genes_use <- VariableFeatures(seurat_obj)
  refdata <- GetAssayData(seurat_obj, assay = "Corrected", slot = "data")[genes_use,]
  
  imputation <- TransferData(
    anchorset = transfer_anchors,
    refdata = refdata,
    weight.reduction = signac_obj[["lsi"]],
    dims = 2:13
  )
  
  signac_obj[["RNA_subclust_imputed"]] <- imputation
  
  return(signac_obj)
}

coembed <- function(seurat_obj, signac_obj){
  signac_obj[["Corrected"]] <- signac_obj[["RNA_subclust_imputed"]]
  signac_obj[["RNA_subclust_imputed"]] <- NULL
  
  genes_use <- VariableFeatures(seurat_obj)
  coembed <- merge(x = seurat_obj[genes_use, ], y = signac_obj[genes_use,])
  coembed <- ScaleData(coembed, features = genes_use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = genes_use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30)
  
  coembed$sequence_type <- "RNA"
  coembed$sequence_type[colnames(coembed) %in% colnames(signac_obj)] <- "ATAC"
  
  return(coembed)
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed_data/clustered_annotated_seurat.qs")
signac <- qread("data/processed_data/signac_integrated.qs")

output_dir <- file.path("data/processed_data/subclustering/")

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
#' go through each major cell type
#' use manual annotation to remove small groups of cells that made it into the wrong cluster
#' then cluster again, integrate, and coembed
#' 

## keratinocytes -------------------------
  ### RNA
kera_rna <- seurat[,seurat$annotation_level1 %in% c("Keratinocyte", "Trichocyte")] 

kera_rna <- cluster_seurat(kera_rna)

DimPlot(kera_rna, group.by = c("sample_ID", "patient_group"))
DimPlot(kera_rna, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(kera_rna, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"), raster = FALSE, label = TRUE)

kera_rna <- kera_rna[,!(kera_rna$Corrected_snn_res.0.2 %in% c(10,11))]
kera_rna <- cluster_seurat(kera_rna)

  ### ATAC
kera_atac <- signac[,signac$annotation_level1 %in% c("Keratinocyte", "Trichocyte")] 
DefaultAssay(kera_atac) <- "ATAC"

kera_atac <- cluster_signac(kera_atac)
DimPlot(kera_atac, group.by = c("sample_ID", "patient_group"))
DimPlot(kera_atac, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(kera_atac, group.by = c("ATAC_snn_res.0.6","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = FALSE)
DimPlot(kera_atac, group.by = c("ATAC_snn_res.0.6","ATAC_snn_res.0.8","ATAC_snn_res.1"), raster = FALSE, label = TRUE)

kera_atac <- kera_atac[,!(kera_atac$ATAC_snn_res.0.6 %in% c(10,9,12))]
kera_atac <- cluster_signac(kera_atac)

  ### integrate and label transfer
kera_atac <- integrate_and_label_transfer(kera_rna, kera_atac)
DimPlot(kera_atac, group.by = c("ATAC_snn_res.0.6","cluster_transfer","ms_veryfine_transfer"), raster = FALSE, label = TRUE)

  ### find cicero connections
# kera_atac <- find_coaccessable_networks(kera_atac) # excluded because it is too slow

  ### coembed
kera_coembed <- coembed(kera_rna, kera_atac)
DimPlot(kera_coembed, group.by = c("sequence_type", "patient_group"))
DimPlot(kera_coembed, group.by = c("patient_group"), split.by = "sequence_type")

  ### save
qsave(kera_rna, paste0(output_dir,"kera_rna.qs"))
qsave(kera_atac, paste0(output_dir,"kera_atac.qs"))
qsave(kera_coembed, paste0(output_dir,"kera_coembed.qs"))




## fibroblasts -------------------------
### RNA
fibro_rna <- seurat[,seurat$annotation_level1 %in% c("Fibroblast")] 

fibro_rna <- cluster_seurat(fibro_rna)

DimPlot(fibro_rna, group.by = c("sample_ID", "patient_group"))
DimPlot(fibro_rna, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(fibro_rna, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"), raster = FALSE, label = TRUE)

fibro_rna <- fibro_rna[,!(fibro_rna$Corrected_snn_res.0.2 %in% c(8))]
fibro_rna <- cluster_seurat(fibro_rna)

### ATAC
fibro_atac <- signac[,signac$annotation_level1 %in% c("Fibroblast")] 
DefaultAssay(fibro_atac) <- "ATAC"

fibro_atac <- cluster_signac(fibro_atac)
DimPlot(fibro_atac, group.by = c("sample_ID", "patient_group"))
DimPlot(fibro_atac, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(fibro_atac, group.by = c("ATAC_snn_res.0.6","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = FALSE)
DimPlot(fibro_atac, group.by = c("ATAC_snn_res.0.6","ATAC_snn_res.0.8","ATAC_snn_res.1"), raster = FALSE, label = TRUE)

# fibro_atac <- fibro_atac[,!(fibro_atac$ATAC_snn_res.0.6 %in% c(10,9,12))]
# fibro_atac <- cluster_signac(fibro_atac)

### integrate and label transfer
fibro_atac <- integrate_and_label_transfer(fibro_rna, fibro_atac)
DimPlot(fibro_atac, group.by = c("ATAC_snn_res.0.6","cluster_transfer","ms_veryfine_transfer"), raster = FALSE, label = TRUE)

### find cicero connections
# fibro_atac <- find_coaccessable_networks(fibro_atac) # excluded because it is too slow

### coembed
fibro_coembed <- coembed(fibro_rna, fibro_atac)
DimPlot(fibro_coembed, group.by = c("sequence_type", "patient_group"))
DimPlot(fibro_coembed, group.by = c("patient_group"), split.by = "sequence_type")

### save
qsave(fibro_rna, paste0(output_dir,"fibro_rna.qs"))
qsave(fibro_atac, paste0(output_dir,"fibro_atac.qs"))
qsave(fibro_coembed, paste0(output_dir,"fibro_coembed.qs"))






## lymphoid -------------------------
### RNA
lymph_rna <- seurat[,seurat$annotation_level1 %in% c("Lymphoid")] 

lymph_rna <- cluster_seurat(lymph_rna)

DimPlot(lymph_rna, group.by = c("sample_ID", "patient_group"))
DimPlot(lymph_rna, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(lymph_rna, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"), raster = FALSE, label = TRUE)

# lymph_rna <- lymph_rna[,!(lymph_rna$Corrected_snn_res.0.2 %in% c(10,11))]
# lymph_rna <- cluster_seurat(lymph_rna)

### ATAC
lymph_atac <- signac[,signac$annotation_level1 %in% c("Lymphoid")] 
DefaultAssay(lymph_atac) <- "ATAC"

lymph_atac <- cluster_signac(lymph_atac)
DimPlot(lymph_atac, group.by = c("sample_ID", "patient_group"))
DimPlot(lymph_atac, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(lymph_atac, group.by = c("ATAC_snn_res.0.6","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = FALSE)
DimPlot(lymph_atac, group.by = c("ATAC_snn_res.0.6","ATAC_snn_res.0.8","ATAC_snn_res.1"), raster = FALSE, label = TRUE)

lymph_atac <- lymph_atac[,!(lymph_atac$ATAC_snn_res.0.6 %in% c(18))]
lymph_atac <- cluster_signac(lymph_atac)

### integrate and label transfer
lymph_atac <- integrate_and_label_transfer(lymph_rna, lymph_atac)
DimPlot(lymph_atac, group.by = c("ATAC_snn_res.0.6","cluster_transfer","ms_veryfine_transfer"), raster = FALSE, label = TRUE)

### find cicero connections
# lymph_atac <- find_coaccessable_networks(lymph_atac) # excluded because it is too slow

### coembed
lymph_coembed <- coembed(lymph_rna, lymph_atac)
DimPlot(lymph_coembed, group.by = c("sequence_type", "patient_group"))
DimPlot(lymph_coembed, group.by = c("patient_group"), split.by = "sequence_type")

### save
qsave(lymph_rna, paste0(output_dir,"lymph_rna.qs"))
qsave(lymph_atac, paste0(output_dir,"lymph_atac.qs"))
qsave(lymph_coembed, paste0(output_dir,"lymph_coembed.qs"))




## Myeloid -------------------------
### RNA
myeloid_rna <- seurat[,seurat$annotation_level1 %in% c("Myeloid","DC")] 

myeloid_rna <- cluster_seurat(myeloid_rna)

DimPlot(myeloid_rna, group.by = c("sample_ID", "patient_group"))
DimPlot(myeloid_rna, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(myeloid_rna, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"), raster = FALSE, label = TRUE)

myeloid_rna <- myeloid_rna[,!(myeloid_rna$Corrected_snn_res.0.2 %in% c(10,11))]
myeloid_rna <- cluster_seurat(myeloid_rna)

### ATAC
myeloid_atac <- signac[,signac$annotation_level1 %in% c("Myeloid","DC")] 
DefaultAssay(myeloid_atac) <- "ATAC"

myeloid_atac <- cluster_signac(myeloid_atac)
DimPlot(myeloid_atac, group.by = c("sample_ID", "patient_group"))
DimPlot(myeloid_atac, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(myeloid_atac, group.by = c("ATAC_snn_res.0.6","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = FALSE)
DimPlot(myeloid_atac, group.by = c("ATAC_snn_res.0.6","ATAC_snn_res.0.8","ATAC_snn_res.1"), raster = FALSE, label = TRUE)

myeloid_atac <- myeloid_atac[,!(myeloid_atac$ATAC_snn_res.0.6 %in% c(10,9,12))]
myeloid_atac <- cluster_signac(myeloid_atac)

### integrate and label transfer
myeloid_atac <- integrate_and_label_transfer(myeloid_rna, myeloid_atac)
DimPlot(myeloid_atac, group.by = c("ATAC_snn_res.0.6","cluster_transfer","ms_veryfine_transfer"), raster = FALSE, label = TRUE)

### find cicero connections
# myeloid_atac <- find_coaccessable_networks(myeloid_atac) # excluded because it is too slow

### coembed
myeloid_coembed <- coembed(myeloid_rna, myeloid_atac)
DimPlot(myeloid_coembed, group.by = c("sequence_type", "patient_group"))
DimPlot(myeloid_coembed, group.by = c("patient_group"), split.by = "sequence_type")

### save
qsave(myeloid_rna, paste0(output_dir,"myeloid_rna.qs"))
qsave(myeloid_atac, paste0(output_dir,"myeloid_atac.qs"))
qsave(myeloid_coembed, paste0(output_dir,"myeloid_coembed.qs"))




## Endothelial -------------------------
### RNA
endo_rna <- seurat[,seurat$annotation_level1 %in% c("Vascular_endo","Lymphatic_endo")] 

endo_rna <- cluster_seurat(endo_rna)

DimPlot(endo_rna, group.by = c("sample_ID", "patient_group"))
DimPlot(endo_rna, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(endo_rna, group.by = c("Corrected_snn_res.0.2","Corrected_snn_res.0.4"), raster = FALSE, label = TRUE)

endo_rna <- endo_rna[,!(endo_rna$Corrected_snn_res.0.2 %in% c(10,11))]
endo_rna <- cluster_seurat(endo_rna)

### ATAC
endo_atac <- signac[,signac$annotation_level1 %in% c("Vascular_endo","Lymphatic_endo")] 
DefaultAssay(endo_atac) <- "ATAC"

endo_atac <- cluster_signac(endo_atac)
DimPlot(endo_atac, group.by = c("sample_ID", "patient_group"))
DimPlot(endo_atac, group.by = c("annotation_level1","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = TRUE)
DimPlot(endo_atac, group.by = c("ATAC_snn_res.0.6","ms_fine_ct", "ms_veryfine_ct"), raster = FALSE, label = FALSE)
DimPlot(endo_atac, group.by = c("ATAC_snn_res.0.6","ATAC_snn_res.0.8","ATAC_snn_res.1"), raster = FALSE, label = TRUE)

endo_atac <- endo_atac[,!(endo_atac$ATAC_snn_res.0.6 %in% c(10,9,12))]
endo_atac <- cluster_signac(endo_atac)

### integrate and label transfer
endo_atac <- integrate_and_label_transfer(endo_rna, endo_atac)
DimPlot(endo_atac, group.by = c("ATAC_snn_res.0.6","cluster_transfer","ms_veryfine_transfer"), raster = FALSE, label = TRUE)

### find cicero connections
# endo_atac <- find_coaccessable_networks(endo_atac) # excluded because it is too slow

### coembed
endo_coembed <- coembed(endo_rna, endo_atac)
DimPlot(endo_coembed, group.by = c("sequence_type", "patient_group"))
DimPlot(endo_coembed, group.by = c("patient_group"), split.by = "sequence_type")

### save
qsave(endo_rna, paste0(output_dir,"endo_rna.qs"))
qsave(endo_atac, paste0(output_dir,"endo_atac.qs"))
qsave(endo_coembed, paste0(output_dir,"endo_coembed.qs"))
















