# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Load functions and libraries                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("utils/rna_utils.R")
set.seed(43648)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_seurat_file <- "data/processed_data/filtered_seurat.qs"

cellranger_folder_path <- "data/rna_GSE212447/GSE212447_RAW/"
file_names_vec <- list.files(cellranger_folder_path)

seurat <- import_seurat(
  cellranger_folder_path = cellranger_folder_path,
  file_names_vec = file_names_vec,
  import_method = "mtx"
)
