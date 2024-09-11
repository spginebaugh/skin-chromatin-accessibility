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
#                                Add Doublets                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frag_path <- "data/atac_GSE212448/fragments/"
frag_files <- list.files("data/atac_GSE212448/fragments/")
frag_files <- frag_files[!(frag_files %in% grep("tbi",frag_files, value = TRUE))]
frag_file_names <- word(frag_files,1,1,"\\.")
frag_file_names <- word(frag_file_names,1,-2,"_")

repeats <- GRanges("chr6", IRanges(1000,2000))
otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))
toExclude <- suppressWarnings(c(repeats, otherChroms))

amulet_doublets <- mclapply(frag_files, function(x) {
  amulet_out <- amulet(paste0(frag_path,x),
                           regionsToExclude=toExclude,
                           fullInMemory = TRUE,
                           BPPARAM = NULL)
                                                       
  return(amulet_out)
}, mc.cores = getOption("mc.cores", 6))
names(amulet_doublets) <- frag_file_names


qsave(amulet_doublets,"data/processed_data/amulet_doublets.qs")
