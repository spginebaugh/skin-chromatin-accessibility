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
library(mclust)
set.seed(43648)

library(ggrastr)

source("utils/archr_helpers.R")
source("utils/matrix_helpers.R")
source("utils/misc_helpers.R")
source("utils/plotting_config.R")


library(BSgenome.Hsapiens.UCSC.hg38)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_file <- "data/processed_data/archRwithpeaks"
atac_proj <- loadArchRProject("data/processed_data/archRsave2")

data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
pointSize <- 0.25


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 call peaks                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Group Coverage Files that can be used for downstream analysis
atac_proj <- addGroupCoverages(
  ArchRProj=atac_proj, 
  groupBy="annotation_level1", 
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40)
  force=TRUE
)


# Find Path to Macs2 binary
pathToMacs2 <- findMacs2()


# Call Reproducible Peaks w/ Macs2
atac_proj <- addReproduciblePeakSet(
  ArchRProj = atac_proj, 
  groupBy = "annotation_level1", 
  peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

# Add Peak Matrix
atac_proj <- addPeakMatrix(ArchRProj = atac_proj, force = TRUE)


# Save project
saveArchRProject(atac_proj, output_file)
