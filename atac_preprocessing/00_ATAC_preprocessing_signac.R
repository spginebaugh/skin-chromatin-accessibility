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
library(scDblFinder)

library(ggrastr)

library(AnnotationHub)
library(BiocParallel)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frag_path <- "data/atac_GSE212448/fragments/"
frag_files <- list.files("data/atac_GSE212448/fragments/")
frag_files <- frag_files[!(frag_files %in% grep("tbi",frag_files, value = TRUE))]
frag_file_names <- word(frag_files,1,1,"\\.")
frag_file_names <- word(frag_file_names,1,-2,"_")

seurat_list <- list()
cutoff <- 1000

for (i in 1:length(frag_files)){
  total_counts_tmp <- CountFragments(paste0(frag_path,frag_files[i]))
  barcodes_tmp <- total_counts_tmp[total_counts_tmp$frequency_count > cutoff, ]$CB
  frags_tmp <- CreateFragmentObject(path = paste0(frag_path,frag_files[i]), cells = barcodes_tmp)
  peaks_tmp <- CallPeaks(frags_tmp)
  counts_tmp <- FeatureMatrix(fragments = frags_tmp, features = peaks_tmp, cells = barcodes_tmp)
  assay_tmp <- CreateChromatinAssay(counts_tmp, fragments = frags_tmp, min.cells = 10, min.features = 200)
  seurat_tmp <- CreateSeuratObject(assay_tmp, assay = "ATAC", project = frag_file_names[i])
    
  seurat_list[[frag_file_names[i]]] <- seurat_tmp
}
qsave(seurat_list, "data/processed_data/signac_list.qs")


seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.id = frag_file_names
)
colnames(seurat@meta.data)[1] <- "sample_ID" ## since we need to demuxlet some of the samples
seurat$barcode <- colnames(seurat)
qsave(seurat,"data/processed_data/signac_merge.qs")












