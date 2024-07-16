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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  functions                               ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
identify_cells <- function(df, TSS_cutoff=6, nFrags_cutoff=2000, minTSS=5, minFrags=1000, maxG=4){
  # Identify likely cells based on gaussian mixture modelling.
  # Assumes that cells, chromatin debris, and other contaminants are derived from
  # distinct gaussians in the TSS x log10 nFrags space. Fit a mixture model to each sample
  # and retain only cells that are derived from a population with mean TSS and nFrags passing
  # cutoffs
  ####################################################################
  # df = data.frame of a single sample with columns of log10nFrags and TSSEnrichment
  # TSS_cutoff = the TSS cutoff that the mean of a generating gaussian must exceed
  # nFrags_cutoff = the log10nFrags cutoff that the mean of a generating gaussian must exceed
  # minTSS = a hard cutoff of minimum TSS for keeping cells, regardless of their generating gaussian
  # maxG = maximum number of generating gaussians allowed
  
  cellLabel <- "cell"
  notCellLabel <- "not_cell"
  
  if(nFrags_cutoff > 100){
    nFrags_cutoff <- log10(nFrags_cutoff)
    minFrags <- log10(minFrags)
  } 
  
  # Fit model
  set.seed(43648)
  mod <- Mclust(df, G=2:maxG, modelNames="VVV")
  
  # Identify classifications that are likely cells
  means <- mod$parameters$mean
  
  # Identify the gaussian with the maximum TSS cutoff
  idents <- rep(notCellLabel, ncol(means))
  idents[which.max(means["TSSEnrichment",])] <- cellLabel
  
  names(idents) <- 1:ncol(means)
  
  # Now return classifications and uncertainties
  df$classification <- idents[mod$classification]
  df$classification[df$TSSEnrichment < minTSS] <- notCellLabel
  df$classification[df$nFrags < minFrags] <- notCellLabel
  df$cell_uncertainty <- NA
  df$cell_uncertainty[df$classification == cellLabel] <- mod$uncertainty[df$classification == cellLabel]
  return(list(results=df, model=mod))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_file <- "data/processed_data/filtered_archr_proj.qs"

proj <- qread("data/processed_data/unfiltered_archr_proj.qs")

sample_cmap <- readRDS("data/utils_data/sample_cmap.rds")
disease_cmap <- head(cmaps_BOR$stallion, 3)
names(disease_cmap) <- c("AA", "C_SD", "C_PB")

pointSize <- 0.5

# Run classification on all samples
minTSS <- 5
samples <- unique(proj$Sample)
cellData <- getCellColData(proj)
cellResults <- lapply(samples, function(x){
  df <- cellData[cellData$Sample == x,c("nFrags","TSSEnrichment")]
  df$log10nFrags <- log10(df$nFrags)
  df <- df[,c("log10nFrags","TSSEnrichment")]
  identify_cells(df, minTSS=minTSS)
})
names(cellResults) <- samples



# # Plot filtering results
# for(samp in samples){
#   df <- as.data.frame(cellResults[[samp]]$results)
#   cell_df <- df[df$classification == "cell",]
#   non_cell_df <- df[df$classification != "cell",]
#   
#   xlims <- c(log10(500), log10(100000))
#   ylims <- c(0, 18)
#   # QC Fragments by TSS plot w/ filtered cells removed:
#   p <- ggPoint(
#     x = cell_df[,1], 
#     y = cell_df[,2], 
#     size = 1.5,
#     colorDensity = TRUE,
#     continuousSet = "sambaNight",
#     xlabel = "Log10 Unique Fragments",
#     ylabel = "TSS Enrichment",
#     xlim = xlims,
#     ylim = ylims,
#     title = sprintf("%s droplets plotted", nrow(cell_df)),
#     rastr = TRUE
#   )
#   # Add grey dots for non-cells
#   p <- p + geom_point_rast(data=non_cell_df, aes(x=log10nFrags, y=TSSEnrichment), color="light grey", size=0.5)
#   p <- p + geom_hline(yintercept = minTSS, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
#   plotPDF(p, name = paste0(samp,"_EM_model_filtered_cells_TSS-vs-Frags.pdf"), ArchRProj = proj, addDOC = FALSE)
# }

finalCellCalls <- lapply(cellResults, function(x) x$results) %>% do.call(rbind, .)
proj <- addCellColData(proj, data=finalCellCalls$classification, name="cellCall", cells=rownames(finalCellCalls), force=TRUE)
proj <- addCellColData(proj, data=finalCellCalls$cell_uncertainty, name="cellCallUncertainty", cells=rownames(finalCellCalls), force=TRUE)


# Add Demuxlet results to C_SD_POOL sample
demuxResults <- c("manuscript_outputs/C_SD_POOL.best")
proj <- addDemuxletResults(proj, bestFiles=demuxResults, sampleNames="C_SD_POOL")

# Relabel demuxlet samples to match existing sample formatting
demuxConvert <- c(
  "C_SD_01_S15" = "C_SD4",
  "C_SD_06_S16" = "C_SD5",
  "C_SD_08_S17" = "C_SD6",
  "C_SD_10_S18" = "C_SD7"
)
proj$Sample2 <- ifelse(proj$DemuxletBest == "NotClassified", proj$Sample, demuxConvert[proj$DemuxletBest])

# Real cells pass QC filter and for C_SD_POOL are classified singlets
realCells <- getCellNames(proj)[(proj$cellCall == "cell") & (proj$DemuxletClassify %ni% c("AMB", "DBL")) & (proj$Sample2 != "GSE212448_C_SD_POOL")]
subProj <- subsetArchRProject(proj, cells=realCells, 
                              outputDirectory="data/filtered_ATAC_output", dropCells=TRUE, force=TRUE)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                add metadata                              ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta_sample <- read_csv("manuscript_metadata/manuscript_table_S1_sample_info.csv")
meta_cell_atac <- read_csv("manuscript_metadata/manuscript_table_S2_scatac_meta.csv")

meta_merge <- data.frame(Sample_ID = subProj$Sample2)
meta_merge <- left_join(meta_merge, meta_sample[,c("Sample_ID","rounded_age","sex","preservation","disease_status")], by = "Sample_ID")


subProj$preservation <- meta_merge$preservation
subProj$sex <- meta_merge$sex %>% factor()
subProj$age <- meta_merge$rounded_age
subProj$disease_status <- meta_merge$disease_status

subProj$diseaseStatus <- NA
subProj$diseaseStatus <- ifelse(grepl("C_SD", subProj$Sample), "C_SD", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("C_PB", subProj$Sample), "C_PB", subProj$diseaseStatus)
subProj$diseaseStatus <- ifelse(grepl("AA", subProj$Sample), "AA", subProj$diseaseStatus)

# Now, add tile matrix and gene score matrix to ArchR project
subProj <- addTileMatrix(subProj, force=TRUE)
subProj <- addGeneScoreMatrix(subProj, force=TRUE)

# Add Infered Doublet Scores to ArchR project (~5-10 minutes)
subProj <- addDoubletScores(subProj, dimsToUse=1:20, scaleDims=TRUE, LSIMethod=2)

# Filter doublets:
subProj <- filterDoublets(subProj, filterRatio = 1)

saveArchRProject(subProj,outputDirectory = "data/processed_data/archRsave1")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  clustering                              ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(43648)

subProj <- addIterativeLSI(
  ArchRProj = subProj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  sampleCellsPre = 15000,
  varFeatures = 50000, 
  dimsToUse = 1:25,
  force = TRUE
)

# Identify Clusters from Iterative LSI
subProj <- addClusters(
  input = subProj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6,
  force = TRUE
)

set.seed(43648)
subProj <- addUMAP(
  ArchRProj = subProj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.4, 
  metric = "cosine",
  force = TRUE
)


# Relabel clusters so they are sorted by cluster size
subProj <- relabelClusters(subProj)

subProj <- addImputeWeights(subProj)

# Make various cluster plots:
subProj <- visualizeClustering(subProj, pointSize=pointSize, sampleCmap=samp_cmap, diseaseCmap=disease_cmap)









