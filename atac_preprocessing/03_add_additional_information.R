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
library(SeuratWrappers)

library(Signac)
library(GenomicRanges)

library(ggrastr)

library(AnnotationHub)
library(BiocParallel)

library(harmony)


library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

library(monocle3)
library(cicero)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seurat <- qread("data/processed_data/signac_clustered.qs")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Motif Info                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(seurat) <- "ATAC"
seqinfo(seurat) <- seqinfo(Annotation(seurat))

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", all_versions = FALSE)  
)

seurat <- AddMotifs(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC",
  pfm = pfm
)

seurat <- RegionStats(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC"
)

seurat <- RunChromVAR(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC"
)

qsave(seurat, "data/processed_data/signac_additional_info.qs")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               Co-accessible network                     ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' cant get conversion from Seurat to CDS to work properly
#' not worth spending the time to fix

# cds <- as.cell_data_set(x = seurat)
# cds <- cluster_cells(cds)
# cds <- as.CellDataSet(cds)
# cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDimS(cds))
# 
# 
# genome <- Signac::seqlengths(seurat)
# genome_df <- data.frame("chr" = names(genome), "length" = genome)
# conns <- run_cicero(cds, genomic_coords = genome.df, sample_num = 100)
# 
# ccans <- generate_ccans(conns)
# 
# links <- ConnectionsToLinks(conns = conns, ccans = ccans)
# 
# ## add into seurat
# Links(seurat) <- links

qsave(seurat, "data/processed_data/signac_all_additional_info.qs")
