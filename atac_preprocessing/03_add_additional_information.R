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
seurat_in <- seurat
rownames(seurat_in) <- str_replace_all(rownames(seurat_in), pattern = "-","_")
cds <- as.cell_data_set(x = seurat)
cds2 <- cds
cds2@assays@data@listData[["counts"]]@Dimnames[[1]] <- str_replace_all(cds2@assays@data@listData[["counts"]]@Dimnames[[1]], pattern = "-","_")
cds2@assays@data@listData[["logcounts"]]@Dimnames[[1]] <- str_replace_all(cds2@assays@data@listData[["logcounts"]]@Dimnames[[1]], pattern = "-","_")

# cds <- cluster_cells(cds)
cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP)

genome <- Signac::seqlengths(seurat)
genome_df <- data.frame("chr" = names(genome), "length" = genome)

distance_parameters <- estimate_distance_parameter(cicero, genomic_coords = genome_df)
mean_distance_parameter <- mean(unlist(distance_parameters))

cicero_model <- generate_cicero_models(cicero,
                                       distance_parameter = mean_distance_parameter,
                                       genomic_coords = genome_df)


conns <- assemble_connections(cicero_model)

ccans <- generate_ccans(conns)

links <- ConnectionsToLinks(conns = conns, ccans = ccans)

## add into seurat
Links(seurat) <- links

qsave(seurat, "data/processed_data/signac_all_additional_info.qs")
