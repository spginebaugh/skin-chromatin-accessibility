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
library(ArchR)
set.seed(43648)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  functions                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create_count_matrix <- function(fragpath, cutoff = 1000){
#   total_counts <- CountFragments(fragpath)
#   barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
#   
#   frags <- CreateFragmentObject(path = fragpath, cells = barcodes)
#   peaks <- CallPeaks(frags)
#   counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
#   return(counts)
# }
# 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
## only have fragment files
fragpath <- "data/atac_GSE212448/fragments/GSE212448_C_PB1_fragments.tsv.gz"

counts1 <- create_count_matrix(fragpath)

arrow_files <- createArrowFiles(
  inputFiles = a,
  sampleNames = names(a),
  geneAnno = geneAnno,
  genomeAnno = genomeAnno,
  minTSS = 0, # Don't filter at this point
  minFrags = 1000, # Default is 1000.
  addTileMat = FALSE, # Don't add tile or geneScore matrices yet. Will add them after we filter
  addGeneScoreMat = FALSE
)
