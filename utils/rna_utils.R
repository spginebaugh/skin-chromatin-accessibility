# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Seurat)
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
library(SingleCellExperiment)

library(scds)
library(scater)
library(DoubletFinder)
library(scDblFinder)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Functions                                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## import to merged seurat object
import_seurat <- function(cellranger_folder_path, file_names_vec, file_h5_path = NA, import_method = "h5") {
  
  seurat_list <- list()
  
  for (file_index in 1:length(file_names_vec)) {
    if (import_method == "h5") {
      seurat_data <- Read10X_h5(filename = paste0(
        cellranger_folder_path,
        file_names_vec[file_index],
        file_h5_path
      ))
    } else if (import_method == "mtx") {
      seurat_data <- Read10X(data.dir = paste0(
        cellranger_folder_path,
        file_names_vec[file_index]
      ))
    } else {
      stop("import_method must be either 'h5' or 'mtx' ")
    }
    seurat_list[file_index] <- CreateSeuratObject(
      counts = seurat_data,
      project = file_names_vec[file_index]
    )
  }
  
  names(seurat_list) <- file_names_vec
  
  seurat <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.id = file_names_vec
  )
  
  seurat <- JoinLayers(seurat) # merge samples in Seurat V5
  
  colnames(seurat@meta.data)[1] <- "sample"
  seurat$barcodes <- colnames(seurat)
  
  return(seurat)
}


## run scDblFinder
run_scDblFinder <- function(seurat_split){
  for (sample_index in 1:length(seurat_split)) {
    dbl_out <- scDblFinder(
      GetAssayData(seurat_split[[sample_index]], assay = "RNA", slot = "counts"),
      returnType = "table"
    ) %>% as.data.frame()
    dbl_out$barcodes <- rownames(dbl_out)
    dbl_out <- dbl_out[, c("barcodes", "class")]
    colnames(dbl_out) <- c("barcodes", "scDbl_class")
    
    ## add scoring into seurat metadata
    metadata <- seurat_split[[sample_index]]@meta.data
    metadata <- left_join(metadata, dbl_out, by = "barcodes")
    rownames(metadata) <- metadata$barcodes
    seurat_split[[sample_index]]@meta.data <- metadata
  }
  return(seurat_split)
}



## run scds doublet scoring
run_scds <- function(seurat_split){
  for (sample_index in 1:length(seurat_split)) {
    sce <- as.SingleCellExperiment(seurat_split[[sample_index]])
    
    sce <- cxds(sce, retRes = TRUE, verb = TRUE)
    sce <- bcds(sce, retRes = TRUE, verb = TRUE)
    sce <- cxds_bcds_hybrid(sce, verb = TRUE)
    
    seurat_split[[sample_index]]$cxds_scores <- sce$cxds_score
    seurat_split[[sample_index]]$bcds_scores <- sce$bcds_score
    seurat_split[[sample_index]]$hybrid_scores <- sce$hybrid_score
  }
  return(seurat_split)
}
