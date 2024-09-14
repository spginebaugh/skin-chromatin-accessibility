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
library(parallel)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frag_path <- "data/atac_GSE212448/fragments/"
frag_files <- list.files("data/atac_GSE212448/fragments/")
frag_files <- frag_files[!(frag_files %in% grep("tbi",frag_files, value = TRUE))]
frag_file_names <- word(frag_files,1,1,"\\.")
frag_file_names <- word(frag_file_names,1,-2,"_")


cutoff <- 1000
total_count_list <- mclapply(frag_files, function(x) {
  total_counts_tmp <- CountFragments(paste0(frag_path,x))
  return(total_counts_tmp)
}, mc.cores = getOption("mc.cores",16))

barcodes_list <- mclapply(total_count_list, function(x) {
  barcodes_tmp <- x[x$frequency_count > cutoff, ]$CB
  return(barcodes_tmp)
}, mc.cores = getOption("mc.cores",16))

frags_list <- list()
for (i in 1:length(frag_files)){
 frags_list[[i]] <- CreateFragmentObject(path = paste0(frag_path,frag_files[i]), cells = barcodes_list[[i]])
}
  
peaks_list <- mclapply(frags_list, function(x) {
  peaks_tmp <- CallPeaks(x, verbose = FALSE)
  return(peaks_tmp)
}, mc.cores = getOption("mc.cores", 16))

combined_peaks <- GenomicRanges::reduce(x = do.call(c, peaks_list))
peakwidths <- width(combined_peaks)
combined_peaks <- combined_peaks[peakwidths < 10000 & peakwidths > 20]

counts_list <- mclapply(1:length(frag_files), function(x){
  counts_tmp <- FeatureMatrix(fragments = frags_list[[x]], 
                              features = combined_peaks, 
                              cells = barcodes_list[[x]])
  return(counts_tmp)
}, mc.cores = getOption("mc.cores", 16))

seurat_list <- mclapply(1:length(frag_files), function(x){
  assay_tmp <- CreateChromatinAssay(counts_list[[x]], fragments = frags_list[[x]], min.cells = 0, min.features = 200)
  seurat_tmp <- CreateSeuratObject(assay_tmp, assay = "ATAC", project = frag_file_names[x])
  return(seurat_tmp)
}, mc.cores = getOption("mc.cores", 16))
names(seurat_list) <- frag_file_names

qsave(seurat_list, "data/processed_data/signac_list.qs")



seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.id = frag_file_names
)
colnames(seurat@meta.data)[1] <- "sample_ID" ## since we need to demuxlet some of the samples
seurat$barcode <- colnames(seurat)
qsave(seurat,"data/processed_data/signac_merge.qs")












