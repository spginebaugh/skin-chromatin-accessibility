# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Load functions and libraries                        ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("utils/rna_utils.R")
set.seed(43648)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 Import Data                               ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_seurat_file <- "data/processed_data/filtered_seurat.qs"

cellranger_folder_path <- "data/rna_GSE212447/GSE212447_RAW/"
file_names_vec <- list.files(cellranger_folder_path)

seurat <- import_seurat(
  cellranger_folder_path = cellranger_folder_path,
  file_names_vec = file_names_vec,
  import_method = "mtx"
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              Add in Metadata                             ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## import
meta_sample <- read_csv("manuscript_metadata/manuscript_table_S1_sample_info.csv")
meta_doublet <- read_csv("data/processed_data/doublet_metadata.csv")
meta_cell_rna <- read_csv("manuscript_metadata/manuscript_table_S3_scrna_meta.csv")


## organize
seurat$sample_ID <- word(seurat$sample,2,-1,"_")
seurat$barcode_ID <- paste0(seurat$sample_ID, "_", word(seurat$barcode,-1,-1,"_"))

colnames(meta_sample)[1] <- "sample_ID"

colnames(meta_cell_rna)[1] <- "barcode_ID"
meta_cell_rna$barcode_ID <- paste0(meta_cell_rna$Sample, "_", word(meta_cell_rna$barcode_ID,-1,-1,"_"))
meta_cell_rna <- meta_cell_rna[,c(1,9,10,11)]
colnames(meta_cell_rna) <- c(
  "barcode_ID",
  "ms_broad_ct", # manuscript broad celltype
  "ms_fine_ct",
  "ms_veryfine_ct"
  )

## combine
metadata <- seurat@meta.data

metadata <- left_join(
  metadata, 
  meta_sample[,c("sample_ID","rounded_age","sex","preservation","disease_status")],
  by = "sample_ID"
)

metadata <- left_join(
  metadata,
  meta_doublet[,c("barcode", "scDbl_class", "cxds_score", "bcds_score",
                  "hybrid_score", "DF_score", "DF_classification")],
  by = "barcode"
)

metadata <- left_join(
  metadata,
  meta_cell_rna,
  by = "barcode_ID"
)

## add back to seurat object
rownames(metadata) <- metadata$barcode
seurat@meta.data <- metadata


## add in cell quality information
seurat$mitoRatio <- PercentageFeatureSet(object = seurat, pattern = "^MT-") / 100
seurat$riboRatio <- PercentageFeatureSet(object = seurat, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA") / 100
seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)
