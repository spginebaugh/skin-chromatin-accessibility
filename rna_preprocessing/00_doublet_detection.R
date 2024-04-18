# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Load functions and libraries                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("utils/rna_utils.R")
source("utils/rna_doublet_scoring_functions.R")

set.seed(43648)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cellranger_folder_path = "data/rna_GSE212447/GSE212447_RAW/"
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
  seurat_split[[i]] <- NormalizeData(seurat_split[[i]])
  seurat_split[[i]] <- FindVariableFeatures(seurat_split[[i]])
  seurat_split[[i]] <- ScaleData(seurat_split[[i]])
  seurat_split[[i]] <- RunPCA(seurat_split[[i]])
}


