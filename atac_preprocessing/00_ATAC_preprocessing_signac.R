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
  assay_tmp <- CreateChromatinAssay(counts_tmp, fragments = frags_tmp, min.cells = 10, min.features = 200)
  seurat_tmp <- CreateSeuratObject(assay_tmp, assay = "ATAC", project = frag_file_names[i])
    
  seurat_list[[frag_file_names[i]]] <- seurat_tmp
}
qsave(seurat_list, "data/processed_data/signac_list.qs")


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #                                Add Doublets                              ----
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mcparam <-(MulticoreParam(workers = multicoreWorkers()))
# 
# repeats <- GRanges("chr6", IRanges(1000,2000))
# otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))
# toExclude <- suppressWarnings(c(repeats, otherChroms))
# 
# amulet_doublets <- list()
# for (i in 1:length(frag_files)){
#   amulet_doublets[[i]] <- amulet(paste0(frag_path,frag_files[i]), 
#                             regionsToExclude=toExclude,
#                             fullInMemory = TRUE,
#                             BPPARAM = mcparam)
# }
# qsave(amulet_doublets,"data/processed_data/amulet_doublets.qs")


seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.id = frag_file_names
)
colnames(seurat@meta.data)[1] <- "sample_ID" ## since we need to demuxlet some of the samples
seurat$barcode <- colnames(seurat)
qsave(seurat,"data/processed_data/signac_merge.qs")

# run_scDblFinder_atac <- function(seurat_split){
#   for (sample_index in 1:length(seurat_split)) {
#     dbl_out <- scDblFinder(
#       GetAssayData(seurat_split[[sample_index]], assay = "ATAC", slot = "counts"),
#       returnType = "table",
#       aggregateFeatures = TRUE,
#       nfeatures=25
#     ) %>% as.data.frame()
#     dbl_out$barcode <- rownames(dbl_out)
#     dbl_out <- dbl_out[, c("barcode", "class")]
#     colnames(dbl_out) <- c("barcode", "scDbl_class")
#     
#     ## add scoring into seurat metadata
#     metadata <- seurat_split[[sample_index]]@meta.data
#     metadata <- left_join(metadata, dbl_out, by = "barcode")
#     rownames(metadata) <- metadata$barcode
#     seurat_split[[sample_index]]@meta.data <- metadata
#   }
#   return(seurat_split)
# }
# 
# seurat_list <- run_scDblFinder_atac(seurat_list)


##################
ccd <- getCellColData(atac_proj)
ccd[ , "DemuxletClassify"] <- "NotClassified"
ccd[ , "DemuxletBest"] <- "NotClassified"
best <- data.frame(readr::read_tsv("manuscript_outputs/C_SD_POOL.best"))
classification <- stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[,1]
cellNames <- paste0("GSE212448_C_SD_POOL", "#", best$BARCODE)
idx <- which(cellNames %in% rownames(ccd))
ccd[ cellNames[idx], "DemuxletClassify"] <- ifelse(
  stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[idx,1] %in% c("AMB","DBL"),
  stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[idx,1],
  stringr::str_split(best$BEST, pattern = "-", simplify=TRUE)[idx,2]
)
ccd[ cellNames[idx], "DemuxletBest"] <- gsub("AMB-","",gsub("DBL-","",gsub("SNG-","",best$BEST)))[idx]
##################


add_demuxlet_results <- function(signac_obj_meta, best_file, sample_names){
  best_dat <- readr::read_tsv(best_file) %>% data.frame()
  best_class <- stringr::str_split(best_dat$BEST, pattern = "-", simplify = TRUE)[,1]
  idx <- which(best_dat$BARCODE %in% rownames(signac_obj_meta))
  
  signac_obj_meta$DemuxletClassify <- "NotClassified"
  signac_obj_meta$DemuxletBest <- "NotClassified"
  
  signac_obj_meta[rownames(signac_obj_meta)[idx], "DemuxletClassify"] <- ifelse(
    stringr::str_split(best_dat$BEST, pattern = "-", simplify=TRUE)[idx,1] %in% c("AMB","DBL"),
    stringr::str_split(best_dat$BEST, pattern = "-", simplify=TRUE)[idx,1],
    stringr::str_split(best_dat$BEST, pattern = "-", simplify=TRUE)[idx,2]
  )
  signac_obj_meta[rownames(signac_obj_meta)[idx], "DemuxletBest"] <- gsub("AMB-","",gsub("DBL-","",gsub("SNG-","",best_dat$BEST)))[idx]
  return(signac_obj_meta)
}

pool_meta <- add_demuxlet_results(seurat_list[["GSE212448_C_SD_POOL"]]@meta.data, "manuscript_outputs/C_SD_POOL.best", "GSE212448_C_SD_POOL")

demuxConvert <- c(
  "C_SD_01_S15" = "C_SD4",
  "C_SD_06_S16" = "C_SD5",
  "C_SD_08_S17" = "C_SD6",
  "C_SD_10_S18" = "C_SD7"
)

pool_meta$sample <- ifelse(pool_meta$DemuxletBest == "NotClassified", pool_meta$orig.ident, demuxConvert[pool_meta$DemuxletBest])
  
seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.id = frag_file_names
)
colnames(seurat@meta.data)[1] <- "sample"
seurat$barcode <- colnames(seurat)

qsave(seurat, "data/processed_data/signac_merged.qs")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Add Metadata                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

## dont have necessary metadata for this step
# seurat$pct_reads_in_peaks <- seurat$peak_region_fragments / seurat$passed_filters * 100

seurat$blacklist_ratio <- FractionCountsInRegion(
  object = seurat, 
  assay = 'ATAC',
  regions = blacklist_hg38_unified
)

qsave(seurat, "data/processed_data/signac_merged2.qs")


DensityScatter(seurat, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
seurat$nucleosome_group <- ifelse(seurat$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = seurat, group.by = 'nucleosome_group')
VlnPlot(
  object = seurat,
  features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Gene Activity                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene.activities <- GeneActivity(seurat)
seurat[['ATAC_RNA']] <- CreateAssayObject(counts = gene.activities)
seurat <- NormalizeData(
  object = seurat,
  assay = 'ATAC_RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)
DefaultAssay(seurat) <- 'ATAC_RNA'

qsave(seurat, "data/processed_data/signac_clustered.qs")












