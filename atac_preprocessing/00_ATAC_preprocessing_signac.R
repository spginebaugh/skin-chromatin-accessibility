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

library(ggrastr)

library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationHub)

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
  assay_tmp <- CreateChromatinAssay(counts_tmp, fragments = frags_tmp)
  seurat_tmp <- CreateSeuratObject(assay_tmp, assay = "ATAC", project = frag_file_names[i])
    
  seurat_list[[frag_file_names[i]]] <- seurat_tmp
}

seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.id = file_names_vec
)
colnames(seurat@meta.data)[1] <- "sample"
seurat$barcode <- colnames(seurat)

qsave(seurat, "data/processed_data/signac_merged.qs")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               Add Annotation                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(seurat) <- annotations

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   QC                                      ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- NucleosomeSignal(object = seurat)
seurat <- TSSEnrichment(object = seurat)
seurat$pct_reads_in_peaks <- seurat$peak_region_fragments / seurat$passed_filters * 100
seurat$blacklist_ratio <- FractionCountsInRegion(
  object = seurat, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)




DensityScatter(seurat, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
seurat$nucleosome_group <- ifelse(seurat$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = seurat, group.by = 'nucleosome_group')
VlnPlot(
  object = seurat,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Filtering                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- subset(
  x = seurat,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             Normalization                                ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = 'q0')
seurat <- RunSVD(seurat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Clustering                                  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DepthCor(seurat)

seurat <- RunUMAP(object = seurat, reduction = 'lsi', dims = 2:30)
seurat <- FindNeighbors(object = seurat, reduction = 'lsi', dims = 2:30)
seurat <- FindClusters(object = seurat, verbose = FALSE, algorithm = 3)
DimPlot(object = seurat, label = TRUE) + NoLegend()



















