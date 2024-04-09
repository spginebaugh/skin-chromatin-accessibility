# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Load Libraries                              ----
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Functions                                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## run scDblFinder -returns split seurat object
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



## run scds doublet scoring -returns split seurat object
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


## run doubletfinder scoring -returns split seurat object
run_doubletfinder <- function(seurat_split){
  ## pK Identification (no ground-truth)
  gc()
  sweep_res <- list()
  sweep_stats <- list()
  pk_out <- list()
  for (sample_index in 1:length(seurat_split)) {
    sweep_res[[sample_index]] <- paramSweep(seurat_split[[sample_index]], PCs = 1:10, sct = FALSE)
    sweep_stats[[sample_index]] <- summarizeSweep(sweep_res[[sample_index]], GT = FALSE)
    pk_out[[sample_index]] <- find.pK(sweep_stats[[sample_index]])
  }
  
  ## Homotypic Doublet Proportion Estimate
  nexp_poi <- list()
  for (sample_index in 1:length(seurat_split)) {
    nexp_poi[[sample_index]] <- round(0.05 * ncol(seurat_split[[sample_index]])) ## Assuming 5.0% doublet formation rate
  }
  
  ## Run DoubletFinder
  for (sample_index in 1:length(seurat_split)) {
    seurat_split[[sample_index]] <- doubletFinder(
      seurat_split[[sample_index]],
      PCs = 1:10, pN = 0.25,
      pK = pk_out[[sample_index]]$pK[which.max(pk_out[[sample_index]]$BCmetric)] %>% as.character() %>% as.numeric(),
      nExp = nexp_poi[[sample_index]], reuse.pANN = FALSE, sct = FALSE
    )
  }
  return(seurat_split)
}

## add doublet metadata to merged seruat object
add_doublet_metadata <- function(seurat_merged, seurat_split){
  metalist <- lapply(seurat_split, function(x) {
    metadata <- x@meta.data
    return(metadata)
  })
  metalist <- lapply(metalist, function(x) {
    colnames(x)[(ncol(x) - 1):ncol(x)] <- c("DF_score", "DF_classification")
    return(x)
  })
  
  metalist <- do.call("rbind", metalist[1:length(metalist)])
  
  seurat_merged$scDbl_class <- metalist$scDbl_class
  seurat_merged$cxds_scores <- metalist$cxds_score
  seurat_merged$bcds_scores <- metalist$bcds_score
  seurat_merged$hybrid_scores <- metalist$hybrid_score
  seurat_merged$DF_score <- metalist$DF_score
  seurat_merged$DF_classification <- metalist$DF_classification
  
  return(seurat_merged)
}


