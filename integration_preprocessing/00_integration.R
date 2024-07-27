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