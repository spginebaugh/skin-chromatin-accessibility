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

seurat$sample_name <- str_remove_all(seurat$sample, "^GS.*?_")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Organize Metadata                           ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in names(dblt)){
  dblt[[i]]$barcode <- paste0(i, "_", rownames(dblt[[i]]))
}
dblt <- do.call(rbind, dblt)
rownames(dblt) <- dblt$barcode
colnames(dblt)[1:6] <- paste0("amulet_", colnames(dblt)[1:6])

colnames(meta_sample)[1] <- "sample_name"

seurat$metaname_barcode <- word(seurat$barcode,2,-1,"_")
colnames(meta_cell_atac) <- paste0("manuscript_", colnames(meta_cell_atac))
colnames(meta_cell_atac)[1] <- "metaname_barcode"
meta_cell_atac$metaname_barcode <- str_replace_all(meta_cell_atac$metaname_barcode,"#","_")

metadata <- seurat@meta.data
metadata <- left_join(metadata, dblt, by = "barcode")
metadata <- left_join(metadata, meta_sample, by = "sample_name")
metadata <- left_join(metadata, meta_cell_atac, by = "metaname_barcode")

rownames(metadata) <- metadata$barcode
seurat@meta.data <- metadata
