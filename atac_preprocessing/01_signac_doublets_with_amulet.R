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
#                                Add Doublets                              ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options(MulticoreParam=MulticoreParam(workers=20))
mcparam <- MulticoreParam(workers = 20)

frag_path <- "data/atac_GSE212448/fragments/"
frag_files <- list.files("data/atac_GSE212448/fragments/")
frag_files <- frag_files[!(frag_files %in% grep("tbi",frag_files, value = TRUE))]
frag_file_names <- word(frag_files,1,1,"\\.")
frag_file_names <- word(frag_file_names,1,-2,"_")

repeats <- GRanges("chr6", IRanges(1000,2000))
otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))
toExclude <- suppressWarnings(c(repeats, otherChroms))

amulet_doublets <- list()
for (i in 1:length(frag_files)){
  amulet_doublets[[i]] <- amulet(paste0(frag_path,frag_files[i]),
                            regionsToExclude=toExclude,
                            fullInMemory = TRUE,
                            BPPARAM = mcparam)
}
qsave(amulet_doublets,"data/processed_data/amulet_doublets.qs")