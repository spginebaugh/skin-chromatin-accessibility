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

library(ggrastr)

source("utils/archr_helpers.R")
source("utils/matrix_helpers.R")
source("utils/misc_helpers.R")
source("utils/plotting_config.R")

#' The addGeneIntegrationMatrix function is broken
#' The fix is provided in the issue discussion on github
#' https://github.com/GreenleafLab/ArchR/issues/2136#issuecomment-2191190744
source("utils/fixed_addGeneIntegrationMatrix.R")

library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(43648)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_file <- "data/processed_data/integrated_archR"
atac_proj <- loadArchRProject("data/processed_data/archRwithpeaks/")
rna_proj <- qread("data/processed_data/clustered_annotated_seurat.qs")

data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
pointSize <- 0.25
barwidth <- 0.9

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          Integrate RNA and ATACseq                       ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atac_proj <- addImputeWeights(atac_proj)


set.seed(43648)
atac_proj <- add_gene_integration_matrix(
  ArchRProj = atac_proj, 
  force = TRUE,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  sampleCellsATAC = 25000, # Default for both was 10000
  sampleCellsRNA = 25000,
  nGenes = 3000, # Default was 2000
  seRNA = rna_proj, # seurat object
  addToArrow = TRUE, # add gene expression to Arrow Files (Set to false initially)
  groupRNA = "annotation_level1", # used to determine the subgroupings specified in groupList (for constrained integration) Additionally this groupRNA is used for the nameGroup output of this function.
  nameCell = "RNA_paired_cell", #Name of column where cell from scRNA is matched to each cell
  nameGroup = "NamedClust_RNA", #Name of column where group from scRNA is matched to each cell
  nameScore = "predictedScore" #Name of column where prediction score from scRNA
)


atac_proj <- addCoAccessibility(
  ArchRProj = atac_proj,
  reducedDims = "IterativeLSI",
  k = 100 # Default is 100
)

atac_proj <- addPeak2GeneLinks( 
  ArchRProj = atac_proj,
  reducedDims = "IterativeLSI",
  k = 100 # Default is 100
)

## add motif info
atac_proj <- addMotifAnnotations(atac_proj, motifSet="cisbp", name="Motif", force=TRUE)

# Add background peaks
atac_proj <- addBgdPeaks(atac_proj, force = TRUE)

atac_proj <- addDeviationsMatrix(
  ArchRProj = atac_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(atac_proj, output_file)
