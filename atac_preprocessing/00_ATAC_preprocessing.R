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

library(ArchR)
set.seed(43648)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_file <- "data/processed_data/unfiltered_archr_proj.qs"

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


proj <- ArchRProject(
  ArrowFiles = arrow_files, 
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  outputDirectory = "data/unfiltered_ATAC_output"
)

unlink("./*.arrow")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    save                                  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## creating the arrow files & proj takes a long time to run, so I'm going to save here
qsave(proj, output_file)


