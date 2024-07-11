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
  set.seed(1)
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
