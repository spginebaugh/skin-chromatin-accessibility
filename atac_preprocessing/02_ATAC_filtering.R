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

library(AnnotationHub)
library(BiocParallel)

library(harmony)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed_data/signac_merge.qs")

dblt <- qread("data/processed_data/amulet_doublets.qs")

meta_sample <- read_csv("manuscript_metadata/manuscript_table_S1_sample_info.csv")

meta_cell_atac <- read_csv("manuscript_metadata/manuscript_table_S2_scatac_meta.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                Demultiplex                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_demuxlet_results <- function(signac_obj_meta, best_file, sample_names){
  best_dat <- readr::read_tsv(best_file) %>% data.frame()
  best_class <- stringr::str_split(best_dat$BEST, pattern = "-", simplify = TRUE)[,1]
  cell_names <- paste0(sample_names,"_", best_dat$BARCODE)
  idx <- which(cell_names %in% rownames(signac_obj_meta))
  
  signac_obj_meta$DemuxletClassify <- "NotClassified"
  signac_obj_meta$DemuxletBest <- "NotClassified"
  
  signac_obj_meta[cell_names[idx], "DemuxletClassify"] <- ifelse(
    stringr::str_split(best_dat$BEST, pattern = "-", simplify=TRUE)[idx,1] %in% c("AMB","DBL"),
    stringr::str_split(best_dat$BEST, pattern = "-", simplify=TRUE)[idx,1],
    stringr::str_split(best_dat$BEST, pattern = "-", simplify=TRUE)[idx,2]
  )
  signac_obj_meta[cell_names[idx], "DemuxletBest"] <- gsub("AMB-","",gsub("DBL-","",gsub("SNG-","",best_dat$BEST)))[idx]
  return(signac_obj_meta)
}

seurat@meta.data <- add_demuxlet_results(seurat@meta.data, "manuscript_outputs/C_SD_POOL.best", "GSE212448_C_SD_POOL")

demuxConvert <- c(
  "C_SD_01_S15" = "C_SD4",
  "C_SD_06_S16" = "C_SD5",
  "C_SD_08_S17" = "C_SD6",
  "C_SD_10_S18" = "C_SD7"
)

seurat$sample <- ifelse(seurat$DemuxletBest == "NotClassified", seurat$sample_ID, demuxConvert[seurat$DemuxletBest])
seurat <- seurat[,!(seurat$DemuxletClassify %in% c("AMB","DBL")) & !(seurat$sample %in% "GSE212448_C_SD_POOL")]

seurat$sample <- str_remove_all(seurat$sample, "^GS.*?_")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Organize Metadata                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in names(dblt)){
  dblt[[i]]$barcode <- paste0(i, "_", rownames(dblt[[i]]))
}
dblt <- do.call(rbind, dblt)
rownames(dblt) <- dblt$barcode
colnames(dblt)[1:6] <- paste0("amulet_", colnames(dblt)[1:6])

colnames(meta_sample)[1] <- "sample"

seurat$metaname_barcode <- word(seurat$barcode,2,-1,"_")
colnames(meta_cell_atac) <- paste0("manuscript_", colnames(meta_cell_atac))
colnames(meta_cell_atac)[1] <- "metaname_barcode"
meta_cell_atac$metaname_barcode <- str_replace_all(meta_cell_atac$metaname_barcode,"#","_")

metadata <- seurat@meta.data
metadata <- left_join(metadata, dblt, by = "barcode")
metadata <- left_join(metadata, meta_sample, by = "sample")
metadata <- left_join(metadata, meta_cell_atac, by = "metaname_barcode")

rownames(metadata) <- metadata$barcode
seurat@meta.data <- metadata

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

peaks_keep <- seqnames(granges(seurat)) %in% standardChromosomes(granges(seurat))
seurat <- seurat[as.vector(peaks_keep), ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   QC                                      ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- NucleosomeSignal(object = seurat)
seurat <- TSSEnrichment(object = seurat)
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

metadata <- seurat@meta.data
metadata %>%
  ggplot(aes(color = sample, x = nCount_ATAC, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 400) +
  geom_vline(xintercept = 20000)

metadata %>%
  ggplot(aes(color = sample, x = TSS.enrichment, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 3) 

metadata %>%
  ggplot(aes(color = sample, x = blacklist_ratio, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 0.01) 

metadata %>%
  ggplot(aes(color = sample, x = nucleosome_signal, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 3) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Filtering                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' were missing an important file from the cellranger output
#' so we will have to use the manuscript outputs to help us with filtering
seurat <- seurat[,!is.na(seurat$manuscript_Sample)]

seurat <- subset(
  x = seurat,
  subset = nCount_ATAC > 400 &
    nCount_ATAC < 20000 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3 &
    amulet_q.value > 0.01 # remove doublets
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
#' first component has strongly negative correlation,
#' so we exclude from downstream processing

seurat <- RunHarmony(seurat, group.by.vars = "sample", reduction.use = 'lsi', dims.use = 2:30, assay.use = "ATAC", project.dim = FALSE)
seurat <- RunUMAP(object = seurat, reduction = 'harmony', dims = 1:29)
seurat <- FindNeighbors(object = seurat, reduction = 'harmony', dims = 1:29)
seurat <- FindClusters(object = seurat, 
                       algorithm = 3, 
                       resolution = c(0.6,0.8,1))
DimPlot(object = seurat, label = TRUE) + NoLegend()
DimPlot(object = seurat, group.by = c("sample","manuscript_BroadClust"), label = TRUE) + NoLegend()
DimPlot(seurat, group.by = c("ATAC_snn_res.0.6","ATAC_snn_res.0.8","ATAC_snn_res.1"))

Idents(seurat) <- seurat$ATAC_snn_res.0.6
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Gene Activity                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene.activities <- GeneActivity(seurat)
seurat[['ATAC_RNA']] <- CreateAssayObject(counts = gene.activities)
seurat <- NormalizeData(
  object = seurat,
  assay = 'ATAC_RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_ATAC_RNA)
)
DefaultAssay(seurat) <- 'ATAC_RNA'

qsave(seurat, "data/processed_data/signac_clustered.qs")


