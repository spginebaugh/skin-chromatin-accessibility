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













