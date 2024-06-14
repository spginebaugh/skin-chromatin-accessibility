# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Load functions and libraries                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages({
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
})

source("utils/rna_utils.R")
source("utils/rna_doublet_scoring_functions.R")

set.seed(43648)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_metadata_file <- "data/processed_data/doublet_metadata.csv"
cellranger_folder_path <- "data/rna_GSE212447/GSE212447_RAW/"
file_names_vec <- list.files(cellranger_folder_path)

seurat <- import_seurat(
  cellranger_folder_path = cellranger_folder_path,
  file_names_vec = file_names_vec,
  import_method = "mtx"
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Run doublet detection                         ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat_split <- SplitObject(seurat, split.by = "sample")

for (i in 1:length(seurat_split)) {
  seurat_split[[i]] <- NormalizeData(seurat_split[[i]], verbose = FALSE)
  seurat_split[[i]] <- FindVariableFeatures(seurat_split[[i]], verbose = FALSE)
  seurat_split[[i]] <- ScaleData(seurat_split[[i]], verbose = FALSE)
  seurat_split[[i]] <- RunPCA(seurat_split[[i]], verbose = FALSE)
}

## functions pulled in from 'rna_doublet_scoring_functions.R'
seurat_split <- run_scDblFinder(seurat_split)
seurat_split <- run_scds(seurat_split)
seurat_split <- run_doubletfinder(seurat_split)
seurat <- add_doublet_metadata(seurat, seurat_split)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    save                                  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(seurat@meta.data, output_metadata_file, row.names = FALSE)
