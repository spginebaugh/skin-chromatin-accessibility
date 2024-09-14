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
  seurat$barcode <- colnames(seurat)
  
  return(seurat)
}




score_panglao <- function(seurat_obj){
  # returns a matrix of scaled celltype scores for each cluster 
  panglao <- readr::read_tsv("data/processed/PanglaoDB_markers_27_Mar_2020.tsv")
  panglao <- panglao[panglao$species %in% c("Hs", "Mm Hs"), ]
  panglao <- panglao[panglao$`official gene symbol` %in% rownames(seurat_obj), ]
  panglao <- split(panglao$`official gene symbol`, panglao$`cell type`)
  
  seurat_obj <- AddModuleScore(seurat_obj, panglao)
  
  metadata <- seurat_obj@meta.data
  
  metadata_name_start <- (ncol(metadata) - length(panglao) + 1)
  metadata_name_end <- ncol(metadata)
  
  pang_ann <- metadata[,metadata_name_start:metadata_name_end] 
  colnames(pang_ann)<- names(panglao)
  
  pang_ann <- scale(pang_ann)
  agg_pang <- aggregate(pang_ann, list(seurat_obj@active.ident), mean)
  rownames(agg_pang) <- agg_pang[, 1]
  agg_pang <- agg_pang[, -1]
  agg_pang <- t(agg_pang)
  return(agg_pang)
}








