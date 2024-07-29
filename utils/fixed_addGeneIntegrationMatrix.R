####################################################################
# Gene Integration Matrix Methods
####################################################################

#' Add a GeneIntegrationMatrix to ArrowFiles or an ArchRProject
#' 
#' This function, will integrate multiple subsets of scATAC cells with a scRNA experiment, compute matched scRNA profiles and
#' then store this in each samples ArrowFile.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of a matrix in the `ArchRProject` containing gene scores to be used for RNA integration.
#' @param matrixName The name to use for the output matrix containing scRNA-seq integration to be stored in the `ArchRProject`.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' This `reducedDims` will be used in weighting the transfer of data to scRNA to scATAC. See `Seurat::TransferData` for more info.
#' @param seRNA A `SeuratObject` or a scRNA-seq `SummarizedExperiment` (cell x gene) to be integrated with the scATAC-seq data.
#' @param groupATAC A column name in `cellColData` of the `ArchRProj` that will be used to determine the subgroupings specified in `groupList`.
#' This is used to constrain the integration to occur across biologically relevant groups.
#' @param groupRNA A column name in either `colData` (if `SummarizedExperiment`) or `metadata` (if `SeuratObject`) of `seRNA` that 
#' will be used to determine the subgroupings specified in `groupList`. This is used to constrain the integration to occur across biologically relevant groups.
#' Additionally this groupRNA is used for the `nameGroup` output of this function.
#' @param groupList A list of cell groupings for both ATAC-seq and RNA-seq cells to be used for RNA-ATAC integration.
#' This is used to constrain the integration to occur across biologically relevant groups. The format of this should be a list of groups 
#' with subgroups of ATAC and RNA specifying cells to integrate from both platforms. 
#' For example `groupList` <- list(groupA = list(ATAC = cellsATAC_A, RNA = cellsRNA_A), groupB = list(ATAC = cellsATAC_B, RNA = cellsRNA_B))
#' @param sampleCellsATAC An integer describing the number of scATAC-seq cells to be used for integration. 
#' This number will be evenly sampled across the total number of cells in the ArchRProject.
#' @param sampleCellsRNA An integer describing the number of scRNA-seq cells to be used for integration.
#' @param embeddingATAC A `data.frame` of cell embeddings such as a UMAP for scATAC-seq cells to be used for density sampling. The `data.frame` object
#' should have a row for each single cell described in `row.names` and 2 columns, one for each dimension of the embedding.
#' @param embeddingRNA A `data.frame` of cell embeddings such as a UMAP for scRNA-seq cells to be used for density sampling. The `data.frame` object
#' should have a row for each single cell described in `row.names` and 2 columns, one for each dimension of the embedding.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a
#' correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param plotUMAP A boolean determining whether to plot a UMAP for each integration block.
#' @param UMAPParams The list of parameters to pass to the UMAP function if "plotUMAP = TRUE". See the function `umap` in the uwot package.
#' @param nGenes The number of variable genes determined by `Seurat::FindVariableGenes()` to use for integration.
#' @param useImputation A boolean value indicating whether to use imputation for creating the Gene Score Matrix prior to integration.
#' @param reduction The Seurat reduction method to use for integrating modalities. See `Seurat::FindTransferAnchors()` for possible reduction methods.
#' @param addToArrow A boolean value indicating whether to add the log2-normalized transcript counts from the integrated matched RNA to the Arrow files.
#' @param scaleTo Each column in the integrated RNA matrix will be normalized to a column sum designated by `scaleTo` prior to adding to Arrow files.
#' @param genesUse If desired a character vector of gene names to use for integration instead of determined ones from Seurat::variableGenes.
#' @param nameCell A column name to add to `cellColData` for the predicted scRNA-seq cell in the specified `ArchRProject`. This is useful for identifying which cell was closest to the scATAC-seq cell.
#' @param nameGroup A column name to add to `cellColData` for the predicted scRNA-seq group in the specified `ArchRProject`. See `groupRNA` for more details.
#' @param nameScore A column name to add to `cellColData` for the predicted scRNA-seq score in the specified `ArchRProject`. These scores represent
#' the assignment accuracy of the group in the RNA cells. Lower scores represent ambiguous predictions and higher scores represent precise predictions. 
#' @param transferParams Additional params to be passed to `Seurat::TransferData`.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exists in the given `input`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param ... Additional params to be added to `Seurat::FindTransferAnchors`
#' @export
add_gene_integration_matrix <- function(
    ArchRProj = NULL,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = NULL,
    groupATAC = NULL,
    groupRNA = NULL,
    groupList = NULL,
    sampleCellsATAC = 10000,
    sampleCellsRNA = 10000,
    embeddingATAC = NULL,
    embeddingRNA = NULL,
    dimsToUse = 1:30,
    scaleDims = NULL,
    corCutOff = 0.75,
    plotUMAP = TRUE,
    UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE),
    nGenes = 2000,
    useImputation = TRUE,
    reduction = "cca",
    addToArrow = TRUE,
    scaleTo = 10000,
    genesUse = NULL,
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore",
    transferParams = list(),
    threads = getArchRThreads(),
    verbose = TRUE,
    force = FALSE,
    logFile = createLogFile("addGeneIntegrationMatrix"),
    ...
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = matrixName, name = "matrixName", valid = c("character"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment", "Seurat"))
  .validInput(input = groupATAC, name = "groupATAC", valid = c("character", "null"))
  .validInput(input = groupRNA, name = "groupRNA", valid = c("character"))
  .validInput(input = groupList, name = "groupList", valid = c("list", "null"))
  .validInput(input = sampleCellsATAC, name = "sampleCellsATAC", valid = c("integer", "null"))
  .validInput(input = sampleCellsRNA, name = "sampleCellsRNA", valid = c("integer", "null"))
  .validInput(input = embeddingATAC, name = "embeddingATAC", valid = c("data.frame", "null"))
  .validInput(input = embeddingRNA, name = "embeddingRNA", valid = c("data.frame", "null"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = plotUMAP, name = "plotUMAP", valid = c("boolean"))
  .validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))  
  .validInput(input = nGenes, name = "nGenes", valid = c("integer"))
  .validInput(input = useImputation, name = "useImputation", valid = c("boolean"))
  .validInput(input = reduction, name = "reduction", valid = c("character"))
  .validInput(input = addToArrow, name = "addToArrow", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = genesUse, name = "genesUse", valid = c("character", "null"))
  .validInput(input = nameCell, name = "nameCell", valid = c("character"))
  .validInput(input = nameGroup, name = "nameGroup", valid = c("character"))
  .validInput(input = nameScore, name = "nameScore", valid = c("character"))
  .validInput(input = transferParams, name = "transferParams", valid = c("list"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logDiffTime("Running Seurat's Integration Stuart* et al 2019", tstart, verbose = verbose, logFile = logFile)
  
  .requirePackage("Seurat", source = "cran")
  
  .logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "Input-Parameters", logFile=logFile)
  
  if(is.null(groupList)){ #If null use all cells (blocking will still occur)
    groupList <- SimpleList()
    groupList[[1]] <- SimpleList(
      ATAC = ArchRProj$cellNames,
      RNA = colnames(seRNA)
    )
  }
  
  #########################################################################################
  # 1. Check All ATAC is Accounted For!
  #########################################################################################
  .logDiffTime("Checking ATAC Input", tstart, verbose = verbose, logFile = logFile)
  
  if (useMatrix %ni% getAvailableMatrices(ArchRProj)) {
    .logMessage(paste0("Matrix ", useMatrix, " does not exist in the provided ArchRProject. See available matrix names from getAvailableMatrices()!"), logFile = logFile)
    stop("Matrix name provided to useMatrix does not exist in ArchRProject!")
  }
  
  if(!is.null(groupATAC)){
    dfATAC <- getCellColData(ArchRProj = ArchRProj, select = groupATAC, drop = FALSE)
  }
  nCell <- rep(0, length(ArchRProj$cellNames))
  names(nCell) <- ArchRProj$cellNames
  
  groupList <- lapply(seq_along(groupList), function(x){
    
    ATAC <- groupList[[x]]$ATAC
    
    if(!is.null(groupATAC)){
      
      if(any(ATAC %in% dfATAC[,1])){
        idx <- which(ATAC %in% dfATAC[,1])
        ATAC2 <- rownames(dfATAC)[which(dfATAC[,1] %in% ATAC[idx])]
        if(length(idx) == length(ATAC)){
          ATAC <- ATAC2
        }else{
          ATAC <- c(ATAC[-idx], ATAC2)
        }
      }
      
    }
    
    SimpleList(ATAC = ATAC, RNA = groupList[[x]]$RNA)
    
  }) %>% SimpleList
  
  for(i in seq_along(groupList)){
    nCell[groupList[[i]]$ATAC] <- nCell[groupList[[i]]$ATAC] + 1
  }
  
  if(!all(nCell == 1)){
    .logMessage(paste0("Missing ", length(which(nCell == 0)), " cells. Found ", length(which(nCell > 1))," overlapping cells from ArchRProj in groupList! Cannot have overlapping/missing cells in ATAC input, check 'groupList' argument!"), logFile = logFile)
    stop("Missing ", length(which(nCell == 0)), " cells. Found ", length(which(nCell > 1))," overlapping cells from ArchRProj in groupList! Cannot have overlapping/missing cells in ATAC input, check 'groupList' argument!")
  }
  
  #########################################################################################
  # 2. Check All RNA is a Cell Name 
  #########################################################################################
  .logDiffTime("Checking RNA Input", tstart, verbose = verbose, logFile = logFile)
  
  #Set up RNA
  if(inherits(seRNA, "SummarizedExperiment")){
    seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
    if(groupRNA %ni% colnames(colData(seRNA))){
      .logMessage("groupRNA not in colData of seRNA", logFile = logFile)
      stop("groupRNA not in colData of seRNA")
    }
    seuratRNA$Group <- paste0(colData(seRNA)[, groupRNA, drop = TRUE])
    rm(seRNA)
  }else{
    if(groupRNA %ni% colnames(seRNA@meta.data)){
      .logMessage("groupRNA not in meta.data of Seurat Object", logFile = logFile)
      stop("groupRNA not in meta.data of Seurat Object")
    }
    seuratRNA <- seRNA
    seuratRNA$Group <- paste0(seRNA@meta.data[,groupRNA])
    rm(seRNA)
  }
  
  if("RNA" %in% names(seuratRNA@assays)){
    DefaultAssay(seuratRNA) <- "RNA"
  }else{
    stop("'RNA' is not present in Seurat Object's Assays! Please make sure that this assay is present!")
  }
  gc()
  
  if(!is.null(groupRNA)){
    dfRNA <- DataFrame(row.names = colnames(seuratRNA), Group = seuratRNA$Group)
  }
  
  groupList <- lapply(seq_along(groupList), function(x){
    
    RNA <- groupList[[x]]$RNA
    
    if(!is.null(groupRNA)){
      
      if(any(RNA %in% dfRNA[,1])){
        idx <- which(RNA %in% dfRNA[,1])
        RNA2 <- rownames(dfRNA)[which(dfRNA[,1] %in% RNA[idx])]
        if(length(idx) == length(RNA)){
          RNA <- RNA2
        }else{
          RNA <- c(RNA[-idx], RNA2)
        }
      }
      
    }
    
    SimpleList(ATAC = groupList[[x]]$ATAC, RNA = RNA)
    
  }) %>% SimpleList
  
  cellRNA <- unlist(lapply(groupList, function(x) x$RNA))
  if(!all(cellRNA %in% colnames(seuratRNA))){
    .logMessage("Found cells for RNA not in colnames(seRNA)! Please retry your input!", logFile = logFile)
    stop("Found cells for RNA not in colnames(seRNA)! Please retry your input!")
  }
  
  seuratRNA <- seuratRNA[, unique(cellRNA)]
  seuratRNA <- NormalizeData(object = seuratRNA, verbose = FALSE)
  
  #########################################################################################
  # 3. Create Integration Blocks
  #########################################################################################
  
  #Check Gene Names And Seurat RowNames
  geneDF <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  sumOverlap <- sum(unique(geneDF$name) %in% unique(rownames(seuratRNA)))
  if(sumOverlap < 5){
    stop("Error not enough overlaps (",sumOverlap,") between gene names from gene scores (ArchR) and rna matrix (seRNA)!")
  }
  .logDiffTime(paste0("Found ", sumOverlap, " overlapping gene names from gene scores and rna matrix!"), tstart, verbose = TRUE, logFile = logFile)
  
  .logDiffTime("Creating Integration Blocks", tstart, verbose = verbose, logFile = logFile)
  
  blockList <- SimpleList()
  
  for(i in seq_along(groupList)){
    
    gLi <- groupList[[i]]
    
    #######################################
    # ATAC
    #######################################
    
    if(length(gLi$ATAC) > sampleCellsATAC){
      
      if(!is.null(embeddingATAC)){
        probATAC <- .getDensity(embeddingATAC[gLi$ATAC,1], embeddingATAC[gLi$ATAC,2])$density
        probATAC <- probATAC / max(probATAC)
        cellsATAC <- gLi$ATAC[order(probATAC, decreasing = TRUE)]
      }else{
        cellsATAC <- sample(gLi$ATAC, length(gLi$ATAC))
      }
      
      cutoffs <- lapply(seq_len(1000), function(x) length(gLi$ATAC) / x) %>% unlist
      blockSize <- ceiling(min(cutoffs[order(abs(cutoffs - sampleCellsATAC))[1]] + 1, length(gLi$ATAC)))
      
      #Density Based Blocking
      nBlocks <- ceiling(length(gLi$ATAC) / blockSize)
      
      blocks <- lapply(seq_len(nBlocks), function(x){
        cellsATAC[seq(x, length(cellsATAC), nBlocks)]
      }) %>% SimpleList
      
    }else{
      
      blocks <- list(gLi$ATAC)
    }
    
    #######################################
    # RNA
    #######################################
    
    if(!is.null(embeddingRNA)){
      probRNA <- .getDensity(embeddingRNA[gLi$RNA,1], embeddingRNA[gLi$RNA,2])$density
      probRNA <- probRNA / max(probRNA)
    }else{
      probRNA <- rep(1, length(gLi$RNA))
    }
    
    blockListi <- lapply(seq_along(blocks), function(x){
      
      SimpleList(
        ATAC = blocks[[x]],
        RNA = sample(x = gLi$RNA, size = min(sampleCellsRNA, length(gLi$RNA)) , prob = probRNA)
      )
      
    }) %>% SimpleList
    
    blockList <- c(blockList, blockListi)
    
  }
  rm(groupList)
  
  #########################################################################################
  # 4. Begin Integration
  #########################################################################################
  .logDiffTime("Prepping Interation Data", tstart, verbose = verbose, logFile = logFile)
  
  #Clean Project For Parallel
  subProj <- ArchRProj
  subProj@imputeWeights <- SimpleList()
  
  #Gene Score Info
  geneDF <- .getFeatureDF(getArrowFiles(subProj), useMatrix)
  geneDF <- geneDF[geneDF$name %in% rownames(seuratRNA), , drop = FALSE]
  
  #Re-Index RNA
  splitGeneDF <- S4Vectors::split(geneDF, geneDF$seqnames)
  featureDF <- lapply(splitGeneDF, function(x){
    x$idx <- seq_len(nrow(x))
    return(x)
  }) %>% Reduce("rbind", .)
  dfParams <- data.frame(
    reduction = reduction
  )
  allChr <- unique(featureDF$seqnames)
  
  #Temp File Prefix
  tmpFile <- .tempfile()
  o <- suppressWarnings(file.remove(paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")))
  
  if(threads > 1){
    h5disableFileLocking()
  }
  
  rD <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  
  #Create Output Directory
  outDir1 <- getOutputDirectory(ArchRProj)
  outDir2 <- file.path(outDir1, "RNAIntegration")
  outDir3 <- file.path(outDir2, matrixName)
  dir.create(outDir1, showWarnings = FALSE)
  dir.create(outDir2, showWarnings = FALSE)
  dir.create(outDir3, showWarnings = FALSE)
  prevFiles <- list.files(outDir3, full.names = TRUE)
  prevFiles <- .suppressAll(file.remove(prevFiles))
  
  tstart <- Sys.time()
  
  threads2 <- max(ceiling(threads * 0.75), 1) #A Little Less here for now
  
  .logDiffTime(paste0("Computing Integration in ", length(blockList), " Integration Blocks!"), tstart, verbose = verbose, logFile = logFile)
  
  #Integration
  dfAll <- .safelapply(seq_along(blockList), function(i){
    
    prefix <- sprintf("Block (%s of %s) :", i , length(blockList))
    
    .logDiffTime(sprintf("%s Computing Integration", prefix), tstart, verbose = verbose, logFile = logFile)
    blocki <- blockList[[i]]
    
    #Subset ATAC
    subProj@cellColData <- subProj@cellColData[blocki$ATAC, ]
    subProj@sampleColData <- subProj@sampleColData[unique(subProj$Sample),,drop=FALSE]
    
    #Subset RNA
    subRNA <- seuratRNA[, blocki$RNA]
    
    #Subet RNA
    subRNA <- subRNA[rownames(subRNA) %in% geneDF$name, ]
    
    ##############################################################################################
    #1. Create Seurat RNA and Normalize
    ##############################################################################################
    .logDiffTime(sprintf("%s Identifying Variable Genes", prefix), tstart, verbose = verbose, logFile = logFile)
    subRNA <- FindVariableFeatures(object = subRNA, nfeatures = nGenes, verbose = FALSE)
    subRNA <- ScaleData(object = subRNA, verbose = FALSE)
    if(is.null(genesUse)){
      genesUse <- VariableFeatures(object = subRNA)
    }
    
    ##############################################################################################
    #2. Get Gene Score Matrix and Create Seurat ATAC
    ##############################################################################################
    .logDiffTime(sprintf("%s Getting GeneScoreMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
    mat <- .getPartialMatrix(
      getArrowFiles(subProj), 
      featureDF = geneDF[geneDF$name %in% genesUse,], 
      threads = 1,
      cellNames = subProj$cellNames,
      useMatrix = useMatrix,
      verbose = FALSE
    )
    rownames(mat) <- geneDF[geneDF$name %in% genesUse, "name"]
    .logThis(mat, paste0("GeneScoreMat-Block-",i), logFile=logFile)
    
    #Impute Matrix (its already scaled internally in ArrowFiles)
    if(useImputation){
      .logDiffTime(sprintf("%s Imputing GeneScoreMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
      imputeParams <- list()
      imputeParams$ArchRProj <- subProj
      imputeParams$randomSuffix <- TRUE
      imputeParams$reducedDims <- reducedDims
      imputeParams$dimsToUse <- dimsToUse
      imputeParams$scaleDims <- scaleDims
      imputeParams$corCutOff <- corCutOff
      imputeParams$threads <- 1
      imputeParams$logFile <- logFile
      subProj <- suppressMessages(do.call(addImputeWeights, imputeParams))
      mat <- suppressMessages(imputeMatrix(mat = mat, imputeWeights = getImputeWeights(subProj), verbose = FALSE, logFile = logFile))
      o <- suppressWarnings(file.remove(unlist(getImputeWeights(subProj)[[1]]))) #Clean Up Space
      .logThis(mat, paste0("GeneScoreMat-Block-Impute-",i), logFile=logFile)
    }
    
    #Log-Normalize 
    mat <- log(mat + 1) #use natural log
    # fix here ------------------
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    rownames(mat) <- as.character(rownames(mat)) 
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    ######################
    
    
    seuratATAC <- Seurat::CreateSeuratObject(counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE])
    seuratATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
    
    #Clean Memory
    rm(mat)
    
    #Set Default Assay
    DefaultAssay(seuratATAC) <- "GeneScore"
    seuratATAC <- Seurat::ScaleData(seuratATAC, verbose = FALSE)
    
    ##############################################################################################
    #3. Transfer Anchors  
    ############################################################################################## 
    .logDiffTime(sprintf("%s Seurat FindTransferAnchors", prefix), tstart, verbose = verbose, logFile = logFile)
    transferAnchors <- .retryCatch({ #This sometimes can crash in mclapply so we can just add a re-run parameter
      gc()
      Seurat::FindTransferAnchors(
        reference = subRNA, 
        query = seuratATAC, 
        reduction = reduction, 
        features = genesUse,
        verbose = FALSE,
        ...
      )
    }, maxAttempts = 2, logFile = logFile)
    .logThis(paste0(utils::capture.output(transferAnchors),collapse="\n"), paste0("transferAnchors-",i), logFile=logFile)
    
    ##############################################################################################
    #4. Transfer Data
    ##############################################################################################
    rDSub <- rD[colnames(seuratATAC),,drop=FALSE]
    .logThis(rDSub, paste0("rDSub-", i), logFile = logFile)
    transferParams$anchorset <- transferAnchors
    transferParams$weight.reduction <- CreateDimReducObject(
      embeddings = rDSub, 
      key = "LSI_", 
      assay = DefaultAssay(seuratATAC)
    )
    transferParams$verbose <- FALSE
    transferParams$dims <- seq_len(ncol(rDSub))
    
    #Group
    .logDiffTime(sprintf("%s Seurat TransferData Cell Group Labels", prefix), tstart, verbose = verbose, logFile = logFile)
    transferParams$refdata <- subRNA$Group
    rnaLabels <- do.call(Seurat::TransferData, transferParams)
    
    #RNA Names
    .logDiffTime(sprintf("%s Seurat TransferData Cell Names Labels", prefix), tstart, verbose = verbose, logFile = logFile)
    transferParams$refdata <- colnames(subRNA)
    rnaLabels2 <- do.call(Seurat::TransferData, transferParams)[,1]
    
    if(addToArrow){
      .logDiffTime(sprintf("%s Seurat TransferData GeneMatrix", prefix), tstart, verbose = verbose, logFile = logFile)
      transferParams$refdata <- GetAssayData(subRNA, assay = "RNA", slot = "data")
      gc()
      matchedRNA <- do.call(Seurat::TransferData, transferParams)
      matchedRNA <- matchedRNA@data
    }
    
    #Match results
    matchDF <- DataFrame(
      cellNames = colnames(seuratATAC), 
      predictionScore = rnaLabels$prediction.score.max,
      predictedGroup = rnaLabels$predicted.id,
      predictedCell = rnaLabels2
    )
    rownames(matchDF) <- matchDF$cellNames
    
    .logDiffTime(sprintf("%s Saving TransferAnchors Joint CCA", prefix), tstart, verbose = verbose, logFile = logFile)
    jointCCA <- DataFrame(transferAnchors@object.list[[1]]@reductions$cca@cell.embeddings)
    jointCCA$Assay <- ifelse(endsWith(rownames(jointCCA), "_reference"), "RNA", "ATAC")
    jointCCA$Group <- NA
    jointCCA$Score <- NA
    jointCCA[paste0(colnames(subRNA), "_reference"), "Group"] <- subRNA$Group
    jointCCA[paste0(matchDF$cellNames, "_query"), "Group"] <- matchDF$predictedGroup
    jointCCA[paste0(matchDF$cellNames, "_query"), "Score"] <- matchDF$predictionScore
    .safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
    
    #Clean Memory
    rm(transferParams, transferAnchors)
    gc()
    
    ##############################################################################################
    #5. Add To Temp Hdf5
    ##############################################################################################
    
    if(addToArrow){
      
      .logDiffTime(sprintf("%s Transferring Paired RNA to Temp File", prefix), tstart, verbose = verbose, logFile = logFile)
      
      #Quickly Write to A Temp Hdf5 File Split By Sample to Then Enable Writing to Each Arrow File
      
      tmpFilei <- paste0(tmpFile, "-IntegrationBlock-", i, ".h5")
      o <- h5createFile(tmpFilei)
      sampleNames <- getCellColData(subProj, "Sample")[matchDF$cellNames, ]
      uniqueSamples <- unique(sampleNames)
      matchedRNA <- .safeSubset( #If Rownames disappeared this will catch that!
        mat = matchedRNA, 
        subsetRows = paste0(featureDF$name), 
        subsetCols = matchDF$cellNames
      )
      
      for(z in seq_along(uniqueSamples)){
        
        mat <- matchedRNA[, which(sampleNames == uniqueSamples[z]), drop = FALSE]
        Group <- uniqueSamples[z]
        
        o <- tryCatch({h5delete(tmpFilei, paste0(Group))}, error = function(x){})
        o <- h5createGroup(tmpFilei, paste0(Group))
        
        #Convert Columns to Rle
        j <- Rle(findInterval(seq(mat@x)-1, mat@p[-1]) + 1)
        
        #Info
        lengthRle <- length(j@lengths)
        lengthI <- length(mat@i)
        
        #Create Data Set
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/i"), storage.mode = "integer", 
                                          dims = c(lengthI, 1), level = 0))
        
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jLengths"), storage.mode = "integer", 
                                          dims = c(lengthRle, 1), level = 0))
        
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jValues"), storage.mode = "integer", 
                                          dims = c(lengthRle, 1), level = 0))
        
        o <- .suppressAll(h5createDataset(tmpFilei, paste0(Group, "/x"), storage.mode = "double", 
                                          dims = c(lengthI, 1), level = 0))
        
        #Write Data Set
        o <- .suppressAll(h5write(obj = mat@i + 1, file = tmpFilei, name = paste0(Group,"/i")))
        o <- .suppressAll(h5write(obj = j@lengths, file = tmpFilei, name = paste0(Group,"/jLengths")))
        o <- .suppressAll(h5write(obj = j@values, file = tmpFilei, name = paste0(Group,"/jValues")))
        o <- .suppressAll(h5write(obj = mat@x, file = tmpFilei, name = paste0(Group, "/x")))
        o <- .suppressAll(h5write(obj = colnames(mat), file = tmpFilei, name = paste0(Group, "/cellNames")))
        #Row Names is always the same
        
      }
      
      rm(matchedRNA, mat, j)
      
    }
    
    .logDiffTime(sprintf("%s Completed Integration", prefix), tstart, verbose = verbose, logFile = logFile)
    
    gc()
    
    matchDF$Block <- Rle(i)
    matchDF
    
  }, threads = threads2) %>% Reduce("rbind", .)
  
  ##############################################################################################
  #5. Plot UMAPs for Co-Embeddings from CCA
  ##############################################################################################
  if(plotUMAP){
    
    for(i in seq_along(blockList)){
      
      o <- tryCatch({
        
        prefix <- sprintf("Block (%s of %s) :", i , length(blockList))
        
        .logDiffTime(sprintf("%s Plotting Joint UMAP", prefix), tstart, verbose = verbose, logFile = logFile)
        
        jointCCA <- readRDS(file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
        
        set.seed(1) # Always do this prior to UMAP
        UMAPParams <- .mergeParams(UMAPParams, list(n_neighbors = 40, min_dist = 0.4, metric="cosine", verbose=FALSE))
        UMAPParams$X <- as.data.frame(jointCCA[, grep("CC_", colnames(jointCCA))])
        UMAPParams$ret_nn <- FALSE
        UMAPParams$ret_model <- FALSE
        UMAPParams$n_threads <- 1
        uwotUmap <- tryCatch({
          do.call(uwot::umap, UMAPParams)
        }, error = function(e){
          errorList <- UMAPParams
          .logError(e, fn = "uwot::umap", info = prefix, errorList = errorList, logFile = logFile)
        })
        
        #Add UMAP and Save Again
        jointCCA$UMAP1 <- uwotUmap[,1]
        jointCCA$UMAP2 <- uwotUmap[,2]
        .safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
        
        p1 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = jointCCA$Assay,
          randomize = TRUE, 
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        p2 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = jointCCA$Group, 
          randomize = TRUE,
          size = 0.2,
          title = paste0(prefix, " colored by scRNA Group"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        pdf(file.path(outDir3, paste0("Save-Block", i,"-JointCCA-UMAP.pdf")), width = 12, height = 6, useDingbats = FALSE)
        ggAlignPlots(p1,p2,type="h")
        dev.off()
        
      }, error = function(e){
        
      })
      
    }
    
  }
  
  ##############################################################################################
  #6. Read sub-matrices and store in ArrowFiles
  ##############################################################################################
  
  if(addToArrow){
    
    .logDiffTime("Transferring Data to ArrowFiles", tstart, verbose = verbose, logFile = logFile)
    
    matrixName <- .isProtectedArray(matrixName)
    
    integrationFiles <- paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")
    
    if(!all(file.exists(integrationFiles))){
      .logMessage("Something went wrong with integration as not all temporary files containing integrated RNA exist!", logFile = logFile)
      stop("Something went wrong with integration as not all temporary files containing integrated RNA exist!")
    }
    
    h5list <- .safelapply(seq_along(integrationFiles), function(x){
      h5ls(integrationFiles[x])
    }, threads = threads)
    
    ArrowFiles <- getArrowFiles(ArchRProj)
    allSamples <- names(ArrowFiles)
    
    o <- .safelapply(seq_along(allSamples), function(y){
      
      sample <- allSamples[y]
      
      prefix <- sprintf("%s (%s of %s)", sample, y, length(ArrowFiles))
      
      .logDiffTime(sprintf("%s Getting GeneIntegrationMatrix From TempFiles!", prefix), tstart, verbose = verbose, logFile = logFile)
      
      sampleIF <- lapply(seq_along(h5list), function(x){
        if(any(h5list[[x]]$group==paste0("/",sample))){
          integrationFiles[x]
        }else{
          NULL
        }
      }) %>% unlist
      
      sampleMat <- lapply(seq_along(sampleIF), function(x){
        
        cellNames <- .h5read(sampleIF[x], paste0(sample, "/cellNames"))
        
        mat <- sparseMatrix(
          i = .h5read(sampleIF[x], paste0(sample, "/i"))[,1], 
          j = as.vector(
            Rle(
              .h5read(sampleIF[x], paste0(sample, "/jValues"))[,1], 
              .h5read(sampleIF[x], paste0(sample, "/jLengths"))[,1]
            )
          ), 
          x = .h5read(sampleIF[x], paste0(sample, "/x"))[,1],
          dims = c(nrow(featureDF), length(cellNames))
        )
        colnames(mat) <- cellNames
        
        mat
        
      }) %>% Reduce("cbind", .)
      
      sampleMat@x <- exp(sampleMat@x) - 1 #Back To Counts
      sampleMat <- .normalizeCols(sampleMat, scaleTo = scaleTo) #Scale to 10,000
      sampleMat <- drop0(sampleMat) # Drop 0's
      rownames(sampleMat) <- paste0(featureDF$name)
      sampleMat <- sampleMat[,ArchRProj$cellNames[BiocGenerics::which(ArchRProj$Sample == sample)], drop = FALSE]
      
      ######################################
      # Initialize SP Mat Group
      ######################################
      o <- .createArrowGroup(ArrowFile = ArrowFiles[sample], group = matrixName, force = force)
      
      o <- .initializeMat(
        ArrowFile = ArrowFiles[sample],
        Group = matrixName,
        Class = "double",
        Units = "NormCounts",
        cellNames = colnames(sampleMat),
        params = dfParams,
        featureDF = featureDF,
        force = force
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictionScore"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictionScore")
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedGroup"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictedGroup")
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedCell"], 
        file = ArrowFiles[sample], 
        name = paste0(matrixName, "/Info/predictedCell")
      )
      
      .logDiffTime(sprintf("%s Adding GeneIntegrationMatrix to ArrowFile!", prefix), tstart, verbose = verbose, logFile = logFile)
      
      for(z in seq_along(allChr)){
        
        chrz <- allChr[z]
        
        .logDiffTime(sprintf("Adding GeneIntegrationMatrix to %s for Chr (%s of %s)!", sample, z, length(allChr)), tstart, verbose = FALSE, logFile = logFile)
        
        idz <- BiocGenerics::which(featureDF$seqnames %bcin% chrz)
        matz <- sampleMat[idz, ,drop=FALSE]
        stopifnot(identical(paste0(featureDF$name[idz]), paste0(rownames(matz))))
        
        #Write sparseMatrix to Arrow File!
        o <- .addMatToArrow(
          mat = matz, 
          ArrowFile = ArrowFiles[sample], 
          Group = paste0(matrixName, "/", chrz), 
          binarize = FALSE,
          addColSums = TRUE,
          addRowSums = TRUE,
          addRowVarsLog2 = TRUE,
          logFile = logFile
        )
        
        #Clean Memory
        rm(matz)
        
        if(z %% 3 == 0 | z == length(allChr)){
          gc()
        }
        
      }
      
      0
      
    }, threads = threads)
    
    o <- suppressWarnings(file.remove(integrationFiles))
    
  }
  
  .logDiffTime("Completed Integration with RNA Matrix", tstart, verbose = verbose, logFile = logFile)
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictedCell,
    name = nameCell,
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictedGroup,
    name = nameGroup,
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj, 
    cells = dfAll$cellNames, 
    data = dfAll$predictionScore,
    name = nameScore,
    force = TRUE
  )
  
  .endLogging(logFile = logFile)
  
  return(ArchRProj)
  
}







##########################################################################################
# Validation Methods
##########################################################################################

.validInput <- function(input = NULL, name = NULL, valid = NULL){
  
  valid <- unique(valid)
  
  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }
  
  if(!is.character(name)){
    stop("name must be a character!")
  }
  
  if("null" %in% tolower(valid)){
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }
  
  av <- FALSE
  
  for(i in seq_along(valid)){
    
    vi <- valid[i]
    
    if(vi == "integer" | vi == "wholenumber"){
      
      if(all(is.numeric(input))){
        #https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
        cv <- min(abs(c(input%%1, input%%1-1)), na.rm = TRUE) < .Machine$double.eps^0.5
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "null"){
      
      cv <- is.null(input)
      
    }else if(vi == "bool" | vi == "boolean" | vi == "logical"){
      
      cv <- is.logical(input)
      
    }else if(vi == "numeric"){
      
      cv <- is.numeric(input)
      
    }else if(vi == "vector"){
      
      cv <- is.vector(input)
      
    }else if(vi == "matrix"){
      
      cv <- is.matrix(input)
      
    }else if(vi == "sparsematrix"){
      
      cv <- is(input, "dgCMatrix")
      
    }else if(vi == "character"){
      
      cv <- is.character(input)
      
    }else if(vi == "factor"){
      
      cv <- is.factor(input)
      
    }else if(vi == "rlecharacter"){
      
      cv1 <- is(input, "Rle")
      if(cv1){
        cv <- is(input@values, "factor") || is(input@values, "character")
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "palette"){
      
      cv <- all(.isColor(input))
      
    }else if(vi == "timestamp"){
      
      cv <- is(input, "POSIXct")
      
    }else if(vi == "dataframe" | vi == "data.frame" | vi == "df"){
      
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
      
    }else if(vi == "fileexists"){
      
      cv <- all(file.exists(input))
      
    }else if(vi == "direxists"){
      
      cv <- all(dir.exists(input))
      
    }else if(vi == "granges" | vi == "gr"){
      
      cv <- is(input, "GRanges")
      
    }else if(vi == "grangeslist" | vi == "grlist"){
      
      cv <- .isGRList(input)
      
    }else if(vi == "list" | vi == "simplelist"){
      
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
      
    }else if(vi == "bsgenome"){
      
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text=input))
      }, error = function(e){
        FALSE
      })
      cv <- any(cv1, cv2)
      
    }else if(vi == "se" | vi == "summarizedexperiment"){
      
      cv <- is(input, "SummarizedExperiment")
      
    }else if(vi == "seurat" | vi == "seuratobject"){
      
      cv <- is(input, "Seurat")
      
    }else if(vi == "txdb"){
      
      cv <- is(input, "TxDb")
      
    }else if(vi == "orgdb"){
      
      cv <- is(input, "OrgDb")
      
    }else if(vi == "bsgenome"){
      
      cv <- is(input, "BSgenome")
      
    }else if(vi == "parallelparam"){
      
      cv <- is(input, "BatchtoolsParam")
      
    }else if(vi == "archrproj" | vi == "archrproject"){
      
      cv <- is(input, "ArchRProject")
      ###validObject(input) check this doesnt break anything if we
      ###add it. Useful to make sure all ArrowFiles exist! QQQ
      
    }else{
      
      stop("Validator is not currently supported by ArchR!")
      
    }
    
    if(cv){
      av <- TRUE
      break
    }   
    
  }
  
  if(av){
    
    return(invisible(TRUE))
    
  }else{
    
    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")
    
  }
  
}

#https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
.isWholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

#https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
.isColor <- function(x = NULL){
  unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)))
}

.isDiscrete <- function(x = NULL){
  is.factor(x) || is.character(x) || is.logical(x)
}

.isGRList <- function(x){
  isList <- grepl("list", class(x), ignore.case=TRUE)
  if(!isList){
    FALSE
  }else{
    allGR <- all(unlist(lapply(x, function(x) is(x, "GRanges") )))
    if(allGR){
      TRUE
    }else{
      FALSE
    }
  }
}

#' Get/Validate BSgenome
#' 
#' This function will attempt to get or validate an input as a BSgenome.
#' 
#' @param genome This option must be one of the following: (i) the name of a valid ArchR-supported genome ("hg38", "hg19", or "mm10"),
#' (ii) the name of a `BSgenome` package (for ex. "BSgenome.Hsapiens.UCSC.hg19"), or (iii) a `BSgenome` object.
#' @param masked A boolean describing whether or not to access the masked version of the selected genome. See `BSgenome::getBSgenome()`.
#' @export
validBSgenome <- function(genome = NULL, masked = FALSE){
  
  .validInput(input = genome, name = "genome", valid = c("character", "bsgenome"))
  .validInput(input = masked, name = "masked", valid = c("boolean"))
  
  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    genome <- tryCatch({
      .requirePackage(genome)
      bsg <- eval(parse(text = genome))
      if(inherits(bsg, "BSgenome")){
        return(bsg)
      }else{
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x){
      BSgenome::getBSgenome(genome, masked = masked)
    })  
    return(genome)
  }else{
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }  
}

.validTxDb <- function(TxDb = NULL){
  stopifnot(!is.null(TxDb))
  if(inherits(TxDb, "TxDb")){
    return(TxDb)
  }else if(is.character(TxDb)){
    return(getTxDb(TxDb)) #change
  }else{
    stop("Cannot validate TxDb options are a valid TxDb or character for getTxDb")
  }
}

.validOrgDb <- function(OrgDb = NULL){
  stopifnot(!is.null(OrgDb))
  if(inherits(OrgDb, "OrgDb")){
    return(OrgDb)
  }else if(is.character(OrgDb)){
    return(getOrgDb(OrgDb)) #change
  }else{
    stop("Cannot validate OrgDb options are a valid OrgDb or character for getOrgDb")
  }
}

.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

.validGeneAnnotation <- function(geneAnnotation = NULL){
  
  if(!inherits(geneAnnotation, "SimpleList")){
    if(inherits(geneAnnotation, "list")){
      geneAnnotation <- as(geneAnnotation, "SimpleList")
    }else{
      stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
    }
  }
  if(identical(sort(tolower(names(geneAnnotation))), c("exons", "genes", "tss"))){
    
    gA <- SimpleList()
    gA$genes <- .validGRanges(geneAnnotation[[grep("genes", names(geneAnnotation), ignore.case = TRUE)]])
    gA$exons <- .validGRanges(geneAnnotation[[grep("exons", names(geneAnnotation), ignore.case = TRUE)]])
    gA$TSS <- .validGRanges(geneAnnotation[[grep("TSS", names(geneAnnotation), ignore.case = TRUE)]])
    
  }else{
    stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
  }
  
  gA
  
}

.validGenomeAnnotation <- function(genomeAnnotation = NULL){
  
  if(!inherits(genomeAnnotation, "SimpleList")){
    if(inherits(genomeAnnotation, "list")){
      genomeAnnotation <- as(genomeAnnotation, "SimpleList")
    }else{
      stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    }
  }
  
  if(identical(sort(tolower(names(genomeAnnotation))), c("blacklist", "chromsizes", "genome"))){
    
    gA <- SimpleList()
    gA$blacklist <- .validGRanges(genomeAnnotation[[grep("blacklist", names(genomeAnnotation), ignore.case = TRUE)]])
    if(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]]=="nullGenome"){
      gA$genome <- "nullGenome"
    }else{
      bsg <- validBSgenome(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
      gA$genome <- bsg@pkgname
    }
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])
    
  }else{
    
    stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    
  }
  
  gA
  
}

.validGeneAnnoByGenomeAnno <- function(geneAnnotation, genomeAnnotation){
  
  allSeqs <- unique(paste0(seqnames(genomeAnnotation$chromSizes)))
  
  geneSeqs <- unique(paste0(seqnames(geneAnnotation$genes)))
  if(!all(geneSeqs %in% allSeqs)){
    geneNotIn <- geneSeqs[which(geneSeqs %ni% allSeqs)]
    message("Found Gene Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(geneNotIn, collapse=","))
    geneAnnotation$genes <- .subsetSeqnamesGR(geneAnnotation$genes, names = allSeqs)
  }
  
  exonSeqs <- unique(paste0(seqnames(geneAnnotation$exons)))
  if(!all(exonSeqs %in% allSeqs)){
    exonNotIn <- exonSeqs[which(exonSeqs %ni% allSeqs)]
    message("Found Exon Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(exonNotIn, collapse=","))
    geneAnnotation$exons <- .subsetSeqnamesGR(geneAnnotation$exons, names = allSeqs)
  }
  
  TSSSeqs <- unique(paste0(seqnames(geneAnnotation$TSS)))
  if(!all(TSSSeqs %in% allSeqs)){
    TSSNotIn <- TSSSeqs[which(TSSSeqs %ni% allSeqs)]
    message("Found TSS Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(TSSNotIn, collapse=","))
    geneAnnotation$TSS <- .subsetSeqnamesGR(geneAnnotation$TSS, names = allSeqs)
  }
  
  geneAnnotation
  
}


.validArchRProject <- function(ArchRProj = NULL){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}


####################################
# Log Tools
####################################

#' Set ArchR Logging
#' 
#' This function will set ArchR logging
#'
#' @param useLogs A boolean describing whether to use logging with ArchR.
#' @export
addArchRLogging <- function(useLogs = TRUE){
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")
  message("Setting ArchRLogging = ", useLogs)
  options(ArchR.logging = useLogs)
  return(invisible(0))
}

#' Get ArchR Logging
#' 
#' This function will get ArchR logging
#'
#' @export
getArchRLogging <- function(){
  ArchRLogging <- options()[["ArchR.logging"]]
  if(!is.logical(ArchRLogging)){
    options(ArchR.logging = TRUE)
    return(TRUE)
  }
  ArchRLogging
}

#' Set ArchR Debugging
#' 
#' This function will set ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @param debug A boolean describing whether to use logging with ArchR.
#' @export
addArchRDebugging <- function(debug = FALSE){
  .validInput(input = debug, name = "debug", valid = "boolean")
  message("Setting ArchRDebugging = ", debug)
  options(ArchR.logging = debug)
  return(invisible(0))
}

#' Get ArchR Debugging
#' 
#' This function will get ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @export
getArchRDebugging <- function(){
  ArchRDebugging <- options()[["ArchR.debugging"]]
  if(!is.logical(ArchRDebugging)){
    options(ArchR.debugging = FALSE)
    return(FALSE)
  }
  ArchRDebugging
}

#' Set ArchR Verbosity for Log Messaging
#' 
#' This function will set ArchR logging verbosity.
#'
#' @param verbose A boolean describing whether to printMessages in addition to logging with ArchR.
#' @export
addArchRVerbose <- function(verbose = TRUE){
  .validInput(input = verbose, name = "verbose", valid = "boolean")
  message("Setting addArchRVerbose = ", verbose)
  options(ArchR.verbose = verbose)
  return(invisible(0))
}

#' Set ArchR Verbosity for Log Messaging
#' 
#' This function will get ArchR logging verbosity.
#'
#' @export
getArchRVerbose <- function(){
  ArchRVerbose <- options()[["ArchR.verbose"]]
  if(!is.logical(ArchRVerbose)){
    options(ArchR.verbose = TRUE)
    return(TRUE)
  }
  ArchRVerbose
}

#' Create a Log File for ArchR
#' 
#' This function will create a log file for ArchR functions. If ArchRLogging is not TRUE
#' this function will return NULL.
#'
#' @param name A character string to add a more descriptive name in log file.
#' @param logDir The path to a directory where log files should be written.
#' @export
createLogFile <- function(
    name = NULL,
    logDir = "ArchRLogs",
    useLogs = getArchRLogging()
){
  
  .validInput(input = name, name = "name", valid = "character")
  .validInput(input = logDir, name = "logDir", valid = "character")
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")
  
  if(!useLogs){
    return(NULL)
  }
  dir.create(logDir, showWarnings = FALSE)
  if(is.null(name)){
    logFile <- .tempfile(pattern = "ArchR", fileext = ".log", tmpdir = logDir)
  }else{
    logFile <- .tempfile(pattern = paste0("ArchR-", name), fileext = ".log", tmpdir = logDir)
  }
  logFile
}

.messageDiffTime <- function(...){ #Deprecated
  .logDiffTime(...)
}

.logDiffTime <- function(
    main = "",
    t1 = NULL,
    verbose = TRUE,
    addHeader = FALSE,
    t2 = Sys.time(),
    units = "mins",
    header = "###########",
    tail = "elapsed.",
    precision = 3,
    logFile = NULL,
    useLogs = getArchRLogging()
){
  
  if(verbose){
    
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),precision))
      if(addHeader){
        msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
      }else{
        msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
      }
      if(getArchRVerbose()) message(msg)
    }, error = function(x){
      if(getArchRVerbose()) message("Time Error : ", x)
    })
    
  }
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(!is.null(logFile)){
    if(file.exists(logFile)){
      logStamp <- tryCatch({
        dt <- abs(round(difftime(t2, t1, units = units),precision))
        if(addHeader){
          msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
        }else{
          msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
        }
        cat(paste0(msg,"\n"), file = logFile, append = TRUE)
      }, error = function(x){
        0
      })
    }
  }
  
  return(invisible(0))
  
}

.startLogging <- function(
    logFile = NULL, 
    useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(file.exists(logFile)){
    return(invisible(0))
  }
  
  .getRam <- function(OS = .Platform$OS.type){
    if(grepl("linux", OS, ignore.case = TRUE)){
      ram <- paste0("Linux : ", as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)))
    }else if(grepl("unix", OS, ignore.case = TRUE)){
      ram <- system("/usr/sbin/system_profiler SPHardwareDataType", intern = TRUE)
      ram <- paste0("MAC : ", gsub("Memory:","",gsub(" ","", grep("Memory", ram, value = TRUE))))
    }else{
      ram <- NA
    }
  }
  
  if(getArchRVerbose()) message("ArchR logging to : ", logFile, 
                                "\nIf there is an issue, please report to github with logFile!")
  
  #Begin With
  cat(.ArchRLogo(ascii = "Package", messageLogo = FALSE), file = logFile, append = FALSE) 
  cat("\nLogging With ArchR!\n\n", file = logFile, append = TRUE) 
  cat(paste0("Start Time : ",Sys.time(),"\n\n"), file = logFile, append = TRUE)
  
  #ArchR Info
  cat("------- ArchR Info\n\n", file = logFile, append = TRUE)
  cat(paste0("ArchRThreads = ", getArchRThreads()), file = logFile, append = TRUE)
  tryCatch({
    if(!is.null(getArchRGenome())){
      cat(paste0("\nArchRGenome = ", getArchRGenome()), file = logFile, append = TRUE)
    }
  }, error = function(x){
  })
  cat("\n\n", file = logFile, append = TRUE)
  
  #Add Info
  cat("------- System Info\n\n", file = logFile, append = TRUE)
  cat(paste0("Computer OS = ", .Platform$OS.type), file = logFile, append = TRUE)
  tryCatch({
    cat(paste0("\nTotal Cores = ", detectCores()), file = logFile, append = TRUE)
  }, error = function(x){
  })
  # tryCatch({
  #     cat(paste0("\nTotal RAM = ", .getRam()), file = logFile, append = TRUE)
  # }, error = function(x){
  # })
  cat("\n\n", file = logFile, append = TRUE)
  
  #Session Info
  cat("------- Session Info\n\n", file = logFile, append = TRUE)
  utils::capture.output(sessionInfo(), file = logFile, append = TRUE)
  cat("\n\n------- Log Info\n\n", file = logFile, append = TRUE)
  
  return(invisible(0))
  
}

.logMessage <- function(
    ..., 
    logFile = NULL,
    verbose = FALSE,   
    useLogs = getArchRLogging()
){
  
  if(getArchRVerbose()){
    msg <- utils::capture.output(message(...), type = "message")
    msg <- paste0(msg, collapse = "\n")
  }else{
    msg <- "SuppressedMessaged due to getArchRVerbose() is FALSE!"
  }
  
  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }
  
  if(verbose){
    message(sprintf("%s : %s", Sys.time(), msg))
  }
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)
  
  return(invisible(0))
  
}

.logHeader <- function(
    name = NULL, 
    logFile = NULL,   
    useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
  
  header <- "###########"
  cat(sprintf("\n%s\n%s : %s\n%s\n\n", header, Sys.time(), name, header), file = logFile, append = TRUE)
  
  return(invisible(0))
}

.logStop <- function(
    ..., 
    logFile = NULL,
    useLogs = getArchRLogging()
){
  
  msg <- utils::capture.output(message(...), type = "message")
  msg <- paste0(msg, collapse = "\n")
  
  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }
  
  if(useLogs){
    if(!is.null(logFile)){
      cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)
    }
  }
  
  stop(sprintf("%s\n", msg), call. = FALSE)
  
  return(invisible(0))
  
}

.logError <- function(
    e = NULL,
    fn = NULL,
    info = NULL, 
    errorList = NULL,
    logFile = NULL,   
    useLogs = getArchRLogging(),
    throwError = TRUE,
    debug = getArchRDebugging()
){
  
  header <- "************************************************************"
  
  if(is.null(logFile)){
    useLogs <- FALSE
  }
  
  if(useLogs){
    #To Log File
    cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile), file = logFile, append = TRUE)
    
    utils::capture.output(print(e), file = logFile, append = TRUE)
    
    if(!is.null(errorList)){
      tryCatch({
        #.safeSaveRDS(errorList, "Save-Error.rds")
        .logThis(errorList, name = "errorList", logFile)
      }, error = function(e){
        cat("Error recording errorList", file = logFile, append = TRUE)
      })
    }
    
    cat(sprintf("\n%s\n\n", header), file = logFile, append = TRUE)
  }
  
  #To Console
  cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile))
  
  if(debug){
    if(!is.null(errorList)){
      debugFile <- paste0(gsub("\\.log", "", logFile), "-debug.rds")
      cat(sprintf("\n%s : ArchRDebugging is set to TRUE, DebugFile = %s\n", Sys.time(), debugFile))
      .safeSaveRDS(errorList, debugFile)
    }
  }
  
  print(e)
  
  cat(sprintf("\n%s\n\n", header))
  
  if(throwError) stop("Exiting See Error Above")
  
  return(invisible(0))
  
}


.logThis <- function(
    x = NULL, 
    name = NULL, 
    logFile = NULL, 
    useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(!file.exists(logFile)){
    stop("logFile does not exist! Something may have deleted this file! Exiting...")
  }
  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
  cat(paste0("\n", Sys.time(), " : ", name, ", Class = ", class(x), "\n"), file = logFile, append = TRUE)
  
  if(missing(x)){
    cat("Data is Missing\n\n", file = logFile, append = TRUE)
    return(invisible(0))
  }
  
  if(is.matrix(x)){
    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    
  }else if(is.data.frame(x)){
    
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    
  }else if(is(x, "dgCMatrix") | is(x, "dgeMatrix")){
    
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    cat(paste0(name, ": NonZeroEntries = ", length(x@x), ", EntryRange = [ ", paste0(range(x@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    
  }else if(is(x, "GRanges")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "SummarizedExperiment")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "DataFrame")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "ArchRProj")){
    
    suppressMessages(utils::capture.output(print(proj), file = logFile, append = TRUE))
    
  }else if(is(x, "SimpleList") | is(x, "list")){
    
    for(i in seq_along(x)){
      
      y <- x[[i]]
      
      if(missing(y)){
        next
      }
      
      if(is.matrix(y)){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is.data.frame(y)){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is(y, "dgCMatrix")){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": NonZeroEntries = ", length(y@x), ", EntryRange = [ ", paste0(range(y@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is(y, "SimpleList") | is(y, "list")){
        
        for(j in seq_along(y)){
          
          z <- y[[j]]
          
          if(missing(z)){
            next
          }
          
          if(is.matrix(z)){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is.data.frame(z)){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is(z, "dgCMatrix")){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": NonZeroEntries = ", length(z@x), ", EntrzRange = [ ", paste0(range(z@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is(y, "SimpleList") | is(y, "list")){
            
            #Only print 2x nested lists
            
          }else{
            
            tryCatch({
              cat("\n", file = logFile, append = TRUE)
              cat(paste0(paste0(name,"$", names(y[j])), ": length = ", length(z), "\n"), file = logFile, append = TRUE)
              suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
              cat("\n", file = logFile, append = TRUE)
            }, error = function(q){
            })
            
          }
          
        }
        
      }else{
        
        tryCatch({
          cat("\n", file = logFile, append = TRUE)
          cat(paste0(paste0(name,"$", names(x[i])), ": length = ", length(y), "\n"), file = logFile, append = TRUE)
          suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
          cat("\n", file = logFile, append = TRUE)
        }, error = function(q){
        })
        
      }
      
    }
    
  }else{
    
    tryCatch({
      cat("\n", file = logFile, append = TRUE)
      cat(paste0(name, ": length = ", length(x), "\n"), file = logFile, append = TRUE)
      suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
      cat("\n", file = logFile, append = TRUE)
    }, error = function(q){
    })
    
  }
  
  cat("\n", file = logFile, append = TRUE)
  return(invisible(0))
  
}

.endLogging <- function(
    logFile = NULL, 
    useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  rL <- readLines(logFile)
  o <- tryCatch({
    t1 <- gsub("Start Time : ","", grep("Start Time", rL, ignore.case = TRUE, value = TRUE))
    mn <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "mins"))
    hr <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "hours"))
    cat("\n------- Completed\n\n", file = logFile, append = TRUE)
    cat(paste0("End Time : ",Sys.time(),"\n"), file = logFile, append = TRUE)
    cat(paste0("Elapsed Time Minutes = ", mn), file = logFile, append = TRUE)
    cat(paste0("\nElapsed Time Hours = ", hr), file = logFile, append = TRUE)
    cat("\n\n-------\n\n\n\n", file = logFile, append = TRUE)
    if(getArchRVerbose()) message("ArchR logging successful to : ", logFile)
  }, error = function(x){
  })
  
  # tryCatch({
  #   R.utils::gzip(logFile, paste0(logFile, ".gz"))
  #   message("ArchR logging successful to : ", paste0(logFile, ".gz"))
  # }, error = function(x){
  # })
  
  return(invisible(0))
  
}


##########################################################################################
# Helper Intermediate Methods
##########################################################################################

.mergeParams <- function(paramInput = NULL, paramDefault = NULL){
  for(i in seq_along(paramDefault)){
    if(!(names(paramDefault)[i] %in% names(paramInput))){
      paramInput[[names(paramDefault)[i]]] <- paramDefault[[i]]
    }
  }
  return(paramInput)
}

.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(source) & is.null(installInfo)){
      if(tolower(source) == "cran"){
        installInfo <- paste0('install.packages("',x,'")')
      }else if(tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      }else{
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

##########################################################################################
# Safe saveRDS check
##########################################################################################

.safeSaveRDS <- function(
    object = NULL, 
    file = "", 
    ascii = FALSE, 
    version = NULL, 
    compress = TRUE, 
    refhook = NULL
){
  #Try to save a test data.frame in location
  testDF <- data.frame(a=1,b=2)
  canSave <- suppressWarnings(tryCatch({
    saveRDS(object = testDF, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
    TRUE
  }, error = function(x){
    FALSE
  }))
  if(!canSave){
    dirExists <- dir.exists(dirname(file))
    if(dirExists){
      stop("Cannot saveRDS. File Path : ", file)
    }else{
      stop("Cannot saveRDS because directory does not exist (",dirname(file),"). File Path : ", file)
    }
  }else{
    saveRDS(object = object, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
  }
}

##########################################################################################
# Stat/Summary Methods
##########################################################################################

.computeKNN <- function(
    data = NULL,
    query = NULL,
    k = 50,
    includeSelf = FALSE,
    ...
){
  .validInput(input = data, name = "data", valid = c("dataframe", "matrix"))
  .validInput(input = query, name = "query", valid = c("dataframe", "matrix"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = includeSelf, name = "includeSelf", valid = c("boolean"))
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  .requirePackage("nabor", source = "cran")
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}

.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

.computeROC <- function(labels = NULL, scores = NULL, name="ROC"){
  .calcAUC <- function(TPR = NULL, FPR = NULL){
    # http://blog.revolutionanalytics.com/2016/11/calculating-auc.html
    dFPR <- c(diff(FPR), 0)
    dTPR <- c(diff(TPR), 0)
    out <- sum(TPR * dFPR) + sum(dTPR * dFPR)/2
    return(out)
  }
  labels <- labels[order(scores, decreasing=TRUE)]
  df <- data.frame(
    False_Positive_Rate = cumsum(!labels)/sum(!labels),
    True_Positive_Rate =  cumsum(labels)/sum(labels)
  )
  df$AUC <- round(.calcAUC(df$True_Positive_Rate,df$False_Positive_Rate),3)
  df$name <- name
  return(df)
}

.getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}

.rowScale <- function(mat = NULL, min = NULL, max = NULL){
  if(!is.null(min)){
    rMin <- min
  }else{
    rMin <- matrixStats::rowMins(mat)
  }
  if(!is.null(max)){
    rMax <- max
  }else{
    rMax <- matrixStats::rowMaxs(mat)
  }
  rScale <- rMax - rMin
  matDiff <- mat - rMin
  matScale <- matDiff/rScale
  out <- list(mat=matScale, min=rMin, max=rMax)
  return(out)
}

.quantileCut <- function(x = NULL, lo = 0.025, hi = 0.975, maxIf0 = TRUE){
  q <- quantile(x, probs = c(lo,hi))
  if(q[2] == 0){
    if(maxIf0){
      q[2] <- max(x)
    }
  }
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

.normalizeCols <- function(mat = NULL, colSm = NULL, scaleTo = NULL){
  if(is.null(colSm)){
    colSm <- Matrix::colSums(mat)
  }
  if(!is.null(scaleTo)){
    mat@x <- scaleTo * mat@x / rep.int(colSm, Matrix::diff(mat@p))
  }else{
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
  }
  return(mat)
}

.safeSubset <- function(mat = NULL, subsetRows = NULL, subsetCols = NULL){
  
  if(!is.null(subsetRows)){
    idxNotIn <- which(subsetRows %ni% rownames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(length(idxNotIn), ncol = ncol(mat)))
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows,]
  }
  
  if(!is.null(subsetCols)){
    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(nrow(mat), ncol = length(idxNotIn)))
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[,subsetCols]
  }
  
  mat
  
}

.groupMeans <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

.groupSums <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowSums(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowSums(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

.groupSds <- function(mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gs <- lapply(unique(groups), function(x){
    if (sparse){
      matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
    }else{
      matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gs) <- unique(groups)
  return(gs)
}

.centerRollMean <- function(v = NULL, k = NULL){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}

##########################################################################################
# Miscellaneous Methods
##########################################################################################

.splitEvery <- function(x = NULL, n = NULL){
  #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  if(is.atomic(x)){
    split(x, ceiling(seq_along(x) / n))
  }else{
    split(x, ceiling(seq_len(nrow(x)) / n))
  }
}

.suppressAll <- function(expr = NULL){
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}

.getAssay <- function(se = NULL, assayName = NULL){
  .assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% .assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se),collapse=", ")))
  }
  return(o)
}

.fileExtension <- function (x = NULL){
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

.checkPath <- function(u = NULL, path = NULL, throwError = TRUE){
  if(is.null(u)){
    out <- TRUE
  }
  out <- lapply(u, function(x, error = TRUE){
    if (Sys.which(x) == "") {
      if(!is.null(path) && file.exists(file.path(path, x))){
        o <- TRUE
      }else{
        if(throwError){
          stop(x, " not found in path, please add ", x, " to path!")
        }else{
          o <- FALSE
        }
      }
    }else{
      o <- TRUE
    }
    return(o)
  }) %>% unlist %>% all
  return(out)
}

.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){
  
  #if the directory doesnt already exist and file.exists evaluates to true, then a file exists with that name
  if(!dir.exists(tmpdir)){
    if(file.exists(tmpdir)){
      stop(paste0("Attempted to create temporary directory ", tmpdir," but a file already exists with this name. Please remove this file and try again!"))
    }
  }
  
  dir.create(tmpdir, showWarnings = FALSE)
  
  if(!dir.exists(tmpdir)){
    stop(paste0("Unable to create temporary directory ", tmpdir,". Check file permissions!")) 
  }
  
  if(addDOC){
    doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }
  
  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))
  
}

.ArchRLogo <- function(ascii = "Logo", messageLogo = TRUE){
  Ascii <- list(
    Package = c("
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    "),
    
    #modified from cyu@athena.mit.edu
    Logo = c("
                                                   / |
                                                 /    \\\
            .                                  /      |.
            \\\\\\                              /        |.
              \\\\\\                          /           `|.
                \\\\\\                      /              |.
                  \\\                    /                |\\\
                  \\\\#####\\\           /                  ||
                ==###########>      /                   ||
                 \\\\##==......\\\    /                     ||
            ______ =       =|__ /__                     ||      \\\\\\\
        ,--' ,----`-,__ ___/'  --,-`-===================##========>
       \\\               '        ##_______ _____ ,--,__,=##,__   ///
        ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
        -,____,---'       \\\\####\\\\________________,--\\\\_##,/
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    ")
  )
  
  if(messageLogo){
    message(Ascii[[ascii]])
  }else{
    Ascii[[ascii]]
  }
  
}

##########################################################################################
# Batch Methods
##########################################################################################

.safelapply <- function(..., threads = 1, preschedule = FALSE){
  
  if(tolower(.Platform$OS.type) == "windows"){
    threads <- 1
  }
  
  if(threads > 1){
    .requirePackage("parallel", source = "cran")
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
    
    errorMsg <- list()
    
    for(i in seq_along(o)){ #Make Sure this doesnt explode!
      if(inherits(o[[i]], "try-error")){
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
      }
    }
    
    if(length(errorMsg) != 0){
      
      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)
      
    }
    
  }else{
    
    o <- lapply(...)
    
  }
  
  o
  
}

.batchlapply <- function(args = NULL, sequential = FALSE){
  
  if(is.null(args$tstart)){
    args$tstart <- Sys.time()
  }
  
  #Determine Parallel Backend
  if(inherits(args$parallelParam, "BatchtoolsParam")){
    
    .logStop("Batchtools not yet fully supported please use local parallel threading!", logFile = args$logFile)
    
    .logDiffTime("Batch Execution w/ BatchTools through BiocParallel!", t1 = args$tstart, verbose = TRUE, logFile = args$logFile)
    
    require(BiocParallel)
    
    args$parallelParam <- btParam
    #Unlink registry Directory
    if(dir.exists(args$registryDir)){
      #Clean Up Registry
      unlink(args$registryDir, recursive = TRUE)# Delete registry directory
    }
    
    #Set Up Registry For Runnning
    args$parallelParam$registryargs <- batchtoolsRegistryargs(
      file.dir = args$registryDir,
      work.dir = getwd(),
      packages = character(0L),
      namespaces = character(0L),
      source = character(0L),
      load = character(0L)
    )
    
    #Register
    BPPARAM <- args$parallelParam
    register(BPPARAM)
    
    #Add To Args
    args$BPPARAM <- BPPARAM
    
    if("..." %in% names(args)){
      args["..."] <- NULL
    }
    
    #Run
    args <- args[names(args) %ni% c("threads", "parallelParam", "subThreading")]
    outlist <- do.call(bplapply, args)
    
  }else{
    
    .logDiffTime("Batch Execution w/ safelapply!", t1 = args$tstart, verbose = TRUE, logFile = args$logFile)
    if(sequential){
      args$subThreads <- args$threads
      args$threads <- 1
    }else{
      if(args$threads > length(args$X)){
        args$subThreads <- floor( args$threads / length(args$X) )
        args$threads <- length(args$X)
      }else{
        args$subThreads <- 1
      }
    }
    
    args <- args[names(args) %ni% c("registryDir", "parallelParam", "subThreading")]
    outlist <- do.call(.safelapply, args)
    
  }
  
  return(outlist)
  
}

.retryCatch <- function(expr, ..., maxAttempts = 3, warnAttempts = FALSE, nameFN = "FN", printInfo = NULL, logFile = NULL){
  currentAttempt <- 0
  completed <- FALSE
  while(!completed & currentAttempt <= maxAttempts){
    currentAttempt <- currentAttempt + 1
    if(currentAttempt > 1){
      .logMessage(nameFN, " : Error occured, attempting again (", currentAttempt - 1, " of ", maxAttempts, ")", logFile = logFile)
    }
    ###########################################################
    tryResult <- tryCatch({
      #########################################################
      #Try Catch Statement Here
      if(warnAttempts){
        out <- return(expr)
      }else{
        out <- suppressWarnings(return(expr))
      }
      #########################################################
      list(out = out, completed = TRUE)
    }, error = function(e){
      list(out = e, completed = FALSE)
    }, ...)
    ###########################################################
    completed <- tryResult$completed
  }
  if(!completed){
    .logMessage(nameFN, " : Error occured and could not be resolved after ", maxAttempts, " additional attempts!", logFile = logFile)
    if(!is.null(printInfo)){
      .logMessage("Error occured at ", printInfo, logFile = logFile)
    }
    print(tryResult[[1]])
    stop()
  }
  
  tryResult[[1]]
  
}


##########################################################################################
# Developer Utils
##########################################################################################

.devMode <- function(package = "ArchR"){
  # fn <- unclass(lsf.str(envir = asNamespace(package), all = TRUE))
  # for(i in seq_along(fn)){
  #   tryCatch({
  #     assign(fn[i], paste0(package,':::', fn[i]), envir=globalenv())
  #     #eval(parse(text=paste0(fn[i], paste0('<<-',package,':::'), fn[i])))
  #   }, error = function(x){
  #   })
  # }
  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for(i in seq_along(fn)){
    tryCatch({
      eval(parse(text=paste0(fn[i], '<-ArchR:::', fn[i])))
    }, error = function(x){
    })
  }
}

.convertToPNG <- function(
    ArchRProj = NULL,
    paths = c("QualityControl"),
    recursive = TRUE,
    outDir = "Figures",
    command = "mv"
){
  
  #If error try
  #brew install fontconfig
  
  .requirePackage("pdftools", source = "cran")
  
  if(!is.null(ArchRProj)){
    paths <- c(paths, file.path(getOutputDirectory(ArchRProj), "Plots"))
  }
  
  pdfFiles <- lapply(seq_along(paths), function(i){
    if(recursive){
      dirs <- list.dirs(paths[i], recursive = FALSE, full.names = FALSE)
      if(length(dirs) > 0){
        pdfs <- lapply(seq_along(dirs), function(j){
          list.files(file.path(paths[i], dirs[j]), full.names = TRUE, pattern = "\\.pdf")
        }) %>% unlist
      }else{
        pdfs <- c()
      }
      pdfs <- c(list.files(paths[i], full.names = TRUE, pattern = "\\.pdf"), pdfs)
    }else{
      pdfs <- list.files(paths[i], full.names = TRUE, pattern = "\\.pdf")
    }
    pdfs
  }) %>% unlist
  
  dir.create(outDir, showWarnings = FALSE)
  
  for(i in seq_along(pdfFiles)){
    print(i)
    tryCatch({
      pdf_convert(
        pdfFiles[i], 
        format = "png", 
        pages = NULL, 
        filenames = file.path(outDir, gsub("\\.pdf", "_%d.png",basename(pdfFiles[i]))),
        dpi = 300, 
        opw = "", 
        upw = "", 
        verbose = TRUE
      )
      system(paste0(command, " ", pdfFiles[i], " ", file.path(outDir, basename(pdfFiles[i]))))
    },error=function(x){
      0
    })
  }
  
}

####################################################################
# Hidden Helper Utils for Arrow Files
####################################################################

.validArrow <- function(ArrowFile = NULL){
  o <- h5closeAll()
  if(h5read(ArrowFile,"Class")!="Arrow"){
    warning(
      "This file is not a valid ArrowFile, this most likely is a bug with previous function where the class was not added.\n",
      "To fix your ArrowFiles :\n",
      "\tlapply(getArrowFiles(ArchRProj), function(x) h5write(obj = 'Arrow', file = x, name = 'Class'))",
      "\nThis will be an error in future versions."
    )
  }
  o <- h5closeAll()
  return(ArrowFile)
}

.isProtectedArray <- function(matrixName = NULL, exclude = NULL){
  protectedArrays <- tolower(c("peakmatrix", "tilematrix", "genescorematrix"))
  if(!is.null(exclude)){
    protectedArrays <- protectedArrays[protectedArrays %ni% tolower(exclude)]
  }
  if(tolower(matrixName) %in% protectedArrays){
    stop(sprintf("Error %s cannot be used as this conflicts with another predefined matrix function!", matrixName))
  }
  matrixName
}

.availableArrays <- function(ArrowFiles = NULL, threads = getArchRThreads()){
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  availableArrays <- .safelapply(seq_along(ArrowFiles), function(x){
    groups <- h5ls(ArrowFiles[x]) %>% {.[.$group=="/" & .$otype=="H5I_GROUP","name"]}
    groups <- groups[!grepl("Fragments|Metadata", groups)]
    groups
  }, threads = threads) %>% Reduce("intersect", .)
  o <- h5closeAll()
  return(availableArrays)
}

.availableSeqnames <- function(ArrowFiles = NULL, subGroup = "Fragments", threads = getArchRThreads()){
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  seqList <- .safelapply(seq_along(ArrowFiles), function(x){
    seqnames <- h5ls(ArrowFiles[x]) %>% {.[.$group==paste0("/",subGroup),]$name}
    seqnames <- seqnames[!grepl("Info", seqnames)]
    seqnames
  }, threads = threads)
  if(!all(unlist(lapply(seq_along(seqList), function(x) identical(seqList[[x]],seqList[[1]]))))){
    stop("Not All Seqnames Identical!")
  }
  o <- h5closeAll()
  return(paste0(seqList[[1]]))
}

.availableChr <- function(ArrowFiles = NULL, subGroup = "Fragments"){
  seqnames <- .availableSeqnames(ArrowFiles, subGroup)
  # if(getArchRChrPrefix()){
  #   seqnames <- seqnames[grep("chr", seqnames, ignore.case = TRUE)]
  # }
  if(length(seqnames) == 0){
    stop("No Chr Found in ArrowFiles!")
  }
  return(seqnames)
}

.availableCells <- function(ArrowFile = NULL, subGroup = NULL, passQC = TRUE){
  if(is.null(subGroup)){
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, "Metadata/CellNames")
    if(passQC){
      passQC <- tryCatch({
        h5read(ArrowFile, "Metadata/PassQC")
      }, error = function(x){
        rep(1, length(cellNames))
      })
      cellNames <- cellNames[which(passQC==1)]
    }
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }else{
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, paste0(subGroup, "/Info/CellNames"))
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }
  return(paste0(sampleName,"#",cellNames))
}

.sampleName <- function(ArrowFile = NULL){
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  o <- h5closeAll()
  return(sampleName)
}

.summarizeArrowContent <- function(ArrowFile = NULL){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  h5DF <- h5ls(ArrowFile)
  
  #Re-Organize Content Info
  h5DF <- h5DF[-which(h5DF$group == "/"),]
  groups <- stringr::str_split(h5DF$group, pattern = "/", simplify=TRUE)[,2]
  groupList <- split(h5DF, groups)
  
  #Split Nested Lists
  groupList2 <- lapply(seq_along(groupList), function(x){
    groupDFx <- groupList[[x]]
    groupx <- gsub(paste0("/", names(groupList)[x]),"",groupDFx$group)
    if(all(groupx=="")){
      groupDFx
    }else{
      subDF <- groupDFx[-which(groupx == ""),]
      split(subDF, stringr::str_split(subDF$group, pattern = "/", simplify=TRUE)[,3])
    }
  })
  names(groupList2) <- names(groupList)
  
  
  o <- h5closeAll()
  
  return(groupList2)
  
}

.getMetadata <- function(ArrowFile = NULL){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  arrowMD <- .summarizeArrowContent(ArrowFile)$Metadata
  
  #Which are same dimensions as cell names
  arrowMD <- arrowMD[which(arrowMD$dim == arrowMD$dim[arrowMD$name=="CellNames"]),]
  
  #Load these into a S4 DataFrame
  md <- lapply(seq_len(nrow(arrowMD)), function(x){
    dfx <- DataFrame(h5read(ArrowFile, paste0(arrowMD$group[x],"/",arrowMD$name[x])))
    colnames(dfx) <- arrowMD$name[x]
    dfx
  }) %>% Reduce("cbind", .)
  
  #Correct CellNames
  md$CellNames <- paste0(sampleName,"#",md$CellNames)
  md$Sample <- Rle(sampleName, nrow(md))
  rownames(md) <- md$CellNames
  md <- md[, -which(colnames(md)=="CellNames")]
  md <- md[,order(colnames(md))]
  
  o <- h5closeAll()
  
  return(md)
}

.getFeatureDF <- function(ArrowFiles = NULL, subGroup = "TileMatrix", threads = getArchRThreads()){
  
  threads <- min(threads, length(ArrowFiles))
  
  .helpFeatureDF <- function(ArrowFile = NULL, subGroup = NULL){
    o <- h5closeAll()
    featureDF <- DataFrame(h5read(ArrowFile, paste0(subGroup,"/Info/FeatureDF")))
    featureDF$seqnames <- Rle(as.character(featureDF$seqnames))
    o <- h5closeAll()
    return(featureDF)
  }
  
  fdf <- .helpFeatureDF(ArrowFiles[1], subGroup = subGroup)
  
  if(length(ArrowFiles) > 1){
    ArrowFiles <- ArrowFiles[-1]
    checkIdentical <- .safelapply(seq_along(ArrowFiles), function(x){
      fdfx <- .helpFeatureDF(ArrowFiles[x], subGroup = subGroup)
      identical(fdfx, fdf)
    }, threads = threads) %>% unlist %>% all
    if(!checkIdentical){
      stop("Error not all FeatureDF for asssay is the same!")
    }
  }
  
  #Re-Order for Split Check!
  newOrder <- split(seq_len(nrow(fdf)), fdf$seqnames) %>% {lapply(seq_along(.), function(x) .[[x]])} %>% Reduce("c", .)
  fdf[newOrder,]
  
}

#####################################################################
# Dropping Group From Hdf5 File
#####################################################################
.createArrowGroup <- function(
    ArrowFile = NULL, 
    group = "GeneScoreMatrix", 
    force = FALSE, 
    verbose = FALSE,
    logFile = NULL
){
  
  ArrowInfo <- .summarizeArrowContent(ArrowFile)
  if(group == "Fragments"){ #This shouldnt happen but just in case
    .logMessage(".createArrowGroup : Cannot create Group over Fragments in Arrow!", logFile = logFile)
    stop("Cannot create Group over Fragments in Arrow!")
  }
  
  if(group %in% names(ArrowInfo)){
    #We Should Check How Big it is if it exists
    ArrowGroup <- ArrowInfo[[group]]
    ArrowGroup <- ArrowGroup[names(ArrowGroup) %ni% c("Info")]
    if(length(ArrowGroup) > 0){
      if(!force){
        .logMessage(".createArrowGroup : Arrow Group already exists! Set force = TRUE to continue!", logFile = logFile)
        stop("Arrow Group already exists! Set force = TRUE to continue!")
      }else{
        .logMessage(".createArrowGroup : Arrow Group already exists! Dropping Group from ArrowFile! This will take ~10-30 seconds!", logFile = logFile)
        if(verbose) message("Arrow Group already exists! Dropping Group from ArrowFile! This will take ~10-30 seconds!")
        o <- .dropGroupsFromArrow(ArrowFile = ArrowFile, dropGroups = group, verbose = verbose, logFile = logFile)
        tryCatch({h5createGroup(ArrowFile , group)}, error=function(e){})
        invisible(return(0))
      }
    }
  }else{
    tryCatch({h5createGroup(ArrowFile , group)}, error=function(e){})
    invisible(return(0))
  }
  
}

.dropGroupsFromArrow <- function(
    ArrowFile = NULL, 
    dropGroups = NULL,
    level = 0,
    verbose = FALSE,
    logFile = NULL
){
  
  tstart <- Sys.time()
  
  #Summarize Arrow Content
  ArrowInfo <- .summarizeArrowContent(ArrowFile)
  
  .logMessage(".dropGroupsFromArrow : Initializing Temp ArrowFile", logFile = logFile)
  
  #We need to transfer first
  outArrow <- .tempfile(fileext = ".arrow")
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")
  
  #1. Metadata First
  .logMessage(".dropGroupsFromArrow : Adding Metadata to Temp ArrowFile", logFile = logFile)
  groupName <- "Metadata"
  o <- h5createGroup(outArrow, groupName)
  
  mData <- ArrowInfo[[groupName]]
  
  for(i in seq_len(nrow(mData))){
    h5name <- paste0(groupName, "/", mData$name[i])
    h5write(.h5read(ArrowFile, h5name), file = outArrow, name = h5name)
  }
  
  #2. Other Groups
  .logMessage(".dropGroupsFromArrow : Adding SubGroups to Temp ArrowFile", logFile = logFile)
  
  groupsToTransfer <- names(ArrowInfo)
  groupsToTransfer <- groupsToTransfer[groupsToTransfer %ni% "Metadata"]
  if(!is.null(dropGroups)){
    groupsToTransfer <- groupsToTransfer[tolower(groupsToTransfer) %ni% tolower(dropGroups)]
  }
  
  for(k in seq_along(groupsToTransfer)){
    
    .logDiffTime(paste0("Transferring ", groupsToTransfer[k]), tstart, verbose = verbose, logFile = logFile)
    
    #Create Group
    groupName <- groupsToTransfer[k]
    o <- h5createGroup(outArrow, groupName)
    
    #Sub Data
    mData <- ArrowInfo[[groupName]]
    
    #Get Order Of Sub Groups (Mostly Related to Seqnames)
    seqOrder <- sort(names(mData))
    if(any(grepl("chr", seqOrder))){
      seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
    }
    
    for(j in seq_along(seqOrder)){
      
      if(verbose) message(j, " ", appendLF = FALSE)
      
      #Create Group
      groupJ <- paste0(groupName, "/", seqOrder[j])
      o <- h5createGroup(outArrow, groupJ)
      
      #Sub mData
      mDataj <- mData[[seqOrder[j]]]
      
      #Transfer Components
      for(i in seq_len(nrow(mDataj))){
        h5name <- paste0(groupJ, "/", mDataj$name[i])
        .suppressAll(h5write(.h5read(ArrowFile, h5name), file = outArrow, name = h5name, level = level))
      }
      
    }
    
    gc()
    
    if(verbose) message("")
    
  }
  
  .logMessage(".dropGroupsFromArrow : Move Temp ArrowFile to ArrowFile", logFile = logFile)
  
  rmf <- file.remove(ArrowFile)
  out <- .fileRename(from = outArrow, to = ArrowFile)
  
  .logDiffTime("Completed Dropping of Group(s)", tstart, logFile = logFile, verbose = verbose)
  
  ArrowFile
  
}

.copyArrows <- function(
    inArrows = NULL,
    outArrows = NULL,
    cellsKeep = NULL,
    level = 0,
    verbose = FALSE,
    logFile = NULL,
    threads = 1
){
  
  stopifnot(length(inArrows) == length(outArrows))
  
  unlist(.safelapply(seq_along(inArrows), function(x){
    .copyArrowSingle(
      inArrow = inArrows[x],
      outArrow = outArrows[x],
      cellsKeep = cellsKeep,
      level = level,
      verbose = verbose,
      logFile = logFile
    )
  }, threads = threads))
  
}

.copyArrowSingle <- function(
    inArrow = NULL,
    outArrow = NULL,
    cellsKeep = NULL,
    level = 0,
    verbose = FALSE,
    logFile = NULL
){
  
  tstart <- Sys.time()
  
  #Summarize Arrow Content
  ArrowInfo <- .summarizeArrowContent(inArrow)
  sampleName <- .sampleName(inArrow)
  
  .logMessage(".copyArrow : Initializing Out ArrowFile", logFile = logFile)
  
  #We need to transfer first
  o <- .suppressAll(file.remove(outArrow))
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")
  
  #1. Metadata First
  .logMessage(".copyArrow : Adding Metadata to Out ArrowFile", logFile = logFile)
  groupName <- "Metadata"
  o <- h5createGroup(outArrow, groupName)
  
  mData <- ArrowInfo[[groupName]]
  cellNames <- .h5read(inArrow, "Metadata/CellNames")
  idx <- which(cellNames %in% stringr::str_split(cellsKeep, pattern="#", simplify=TRUE)[,2])
  
  if(length(idx)==0){
    stop("No cells matching in arrow file!")
  }
  
  for(i in seq_len(nrow(mData))){
    h5name <- paste0(groupName, "/", mData$name[i])
    mDatai <- .h5read(inArrow, h5name)
    if(length(mDatai)==length(cellNames)){
      mDatai <- mDatai[idx]
    }
    h5write(mDatai, file = outArrow, name = h5name)
  }
  
  #2. scATAC-Fragments
  .logDiffTime(paste0("Transferring Fragments"), tstart, verbose = verbose, logFile = logFile)
  
  #Create Group
  groupName <- "Fragments"
  o <- h5createGroup(outArrow, groupName)
  
  #Sub Data
  mData <- ArrowInfo[[groupName]]
  
  #Get Order Of Sub Groups (Mostly Related to Seqnames)
  seqOrder <- sort(names(mData))
  if(any(grepl("chr", seqOrder))){
    seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
  }
  
  for(j in seq_along(seqOrder)){
    
    if(verbose) message(j, " ", appendLF = FALSE)
    
    #Create Group
    groupJ <- paste0(groupName, "/", seqOrder[j])
    o <- h5createGroup(outArrow, groupJ)
    
    #Sub mData
    mDataj <- mData[[seqOrder[j]]]
    
    #Read In Fragments
    RGLengths <- .h5read(inArrow, paste0(groupJ, "/RGLengths"))
    RGValues <- .h5read(inArrow, paste0(groupJ, "/RGValues"))
    RGRle <- Rle(paste0(sampleName, "#", RGValues), RGLengths)
    
    #Determine Which to Keep
    idxj <- BiocGenerics::which(RGRle %bcin% cellsKeep)
    
    if(length(idxj) == 0){
      idxj <- 1
    }
    
    #Info
    Ranges <- .h5read(inArrow, paste0(groupJ, "/Ranges"))[idxj, ,drop=FALSE]
    RGRle <- RGRle[idxj]
    RGLengths <- RGRle@lengths
    RGValues <- stringr::str_split(RGRle@values, pattern = "#", simplify = TRUE)[,2]
    
    #Write Barcodes
    o <- .suppressAll(h5write(RGLengths, file = outArrow, name = paste0(groupJ, "/RGLengths"), level = level))
    o <- .suppressAll(h5write(RGValues, file = outArrow, name = paste0(groupJ, "/RGValues"), level = level))
    
    #Write Ranges
    o <- .suppressAll(
      h5write(
        obj = Ranges, 
        file = outArrow, 
        name = paste0(groupJ, "/Ranges"), 
        level = level
      )
    )
    
  }
  
  if(verbose) message("")
  
  #3. Other Matrices
  .logMessage(".copyArrow : Adding SubMatrices to Out ArrowFile", logFile = logFile)
  groupsToTransfer <- names(ArrowInfo)
  groupsToTransfer <- groupsToTransfer[groupsToTransfer %ni% c("Metadata", "Fragments")]
  
  for(k in seq_along(groupsToTransfer)){
    
    .logDiffTime(paste0("Transferring ", groupsToTransfer[k]), tstart, verbose = verbose, logFile = logFile)
    
    #Create Group
    groupName <- groupsToTransfer[k]
    o <- h5createGroup(outArrow, groupName)
    
    #Sub Data
    mData <- ArrowInfo[[groupName]]
    
    #Get Order Of Sub Groups (Mostly Related to Seqnames)
    seqOrder <- sort(names(mData))
    if(any(grepl("chr", seqOrder))){
      seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
    }
    
    cellNames <- paste0(sampleName, "#", .h5read(inArrow, paste0(groupName, "/Info/CellNames")))
    featureDF <- .getFeatureDF(ArrowFile = inArrow, subGroup = groupName)
    seqOrder <- c("Info", seqOrder[!grepl("Info", seqOrder)])
    
    for(j in seq_along(seqOrder)){
      
      if(verbose) message(j, " ", appendLF = FALSE)
      
      #Create Group
      groupJ <- paste0(groupName, "/", seqOrder[j])
      
      if(seqOrder[j] == "Info"){
        
        o <- h5createGroup(outArrow, groupJ)
        
        #Sub mData
        mDataj <- mData[[seqOrder[j]]]
        idxCL <- which(mDataj$dim == mDataj$dim[mDataj$name=="CellNames"])
        idxCL <- idxCL[mDataj$name[idxCL] %ni% "FeatureDF"]
        idxKeep <- which(cellNames %in% cellsKeep)
        
        #Transfer Components
        for(i in seq_len(nrow(mDataj))){
          
          h5name <- paste0(groupJ, "/", mDataj$name[i])
          
          if(i %in% idxCL){
            .suppressAll(h5write(.h5read(inArrow, h5name)[idxKeep], file = outArrow, name = h5name))
          }else{
            .suppressAll(h5write(.h5read(inArrow, h5name), file = outArrow, name = h5name))
          }
          
        }
        
      }else{
        
        #Sub mData
        mDataj <- mData[[seqOrder[j]]]
        addAnalysis <- mDataj[mDataj$name %ni% c("i", "jLengths", "jValues", "x"), "name"]
        
        mat <- .getMatFromArrow(
          ArrowFile = inArrow,
          featureDF = featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqOrder[j]),],
          useMatrix = groupName,
          cellNames = cellNames[cellNames %in% cellsKeep]
        )
        
        o <- .addMatToArrow(
          mat = mat, 
          ArrowFile = outArrow, 
          Group = paste0(groupName, "/", seqOrder[j]), 
          binarize = all(mat@x == 1),
          addColSums = "colSums" %in% addAnalysis,
          addRowSums = "rowSums" %in% addAnalysis,
          addRowMeans = "rowMeans" %in% addAnalysis,
          addRowVars = "rowVars" %in% addAnalysis,
          addRowVarsLog2 = "rowVarsLog2" %in% addAnalysis
        )
        
        rm(mat)
        
      }
      
    }
    
    gc()
    
    if(verbose) message("")
    
  }
  
  .logDiffTime("Completed Copying ArrowFile", tstart, logFile = logFile, verbose = verbose)
  
  outArrow
  
}

##########################################################################################
# Helper Functions for GenomicRanges
##########################################################################################

#' Filters unwanted seqlevels from a Genomic Ranges object or similar object
#'
#' This function allows for removal of manually designated or more broadly undesirable seqlevels from a Genomic Ranges object or similar object
#'
#' @param gr A `GRanges` object or another object containing `seqlevels`.
#' @param remove A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels.
#' If no manual removal is desired, `remove` should be set to `NULL`.
#' @param underscore A boolean value indicating whether to remove all seqlevels whose names contain an underscore (for example "chr11_KI270721v1_random").
#' @param standard A boolean value indicating whether only standard chromosomes should be kept. Standard chromosomes are defined by
#' `GenomeInfoDb::keepStandardChromosomes()`.
#' @param pruningMode The name of the pruning method to use (from`GenomeInfoDb::seqinfo()`) when seqlevels must be removed from a `GRanges` object.
#' When some of the seqlevels to drop from the given `GRanges` object are in use (i.e. have ranges on them), the ranges on these sequences need
#' to be removed before the seqlevels can be dropped. Four pruning modes are currently defined: "error", "coarse", "fine", and "tidy".
#' @export
filterChrGR <- function(
    gr = NULL, 
    remove = NULL, 
    underscore = TRUE, 
    standard = TRUE, 
    pruningMode="coarse"
){
  
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = remove, name = "remove", valid = c("character", "null"))
  .validInput(input = underscore, name = "underscore", valid = c("boolean"))
  .validInput(input = standard, name = "standard", valid = c("boolean"))
  .validInput(input = pruningMode, name = "pruningMode", valid = c("character"))
  
  #first we remove all non standard chromosomes
  if(standard){
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = pruningMode)
  }
  #Then check for underscores or specified remove
  seqNames <- seqlevels(gr)
  chrRemove <- c()
  #first we remove all chr with an underscore
  if(underscore){
    chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
  }
  #next we remove all chr specified in remove
  if(!is.null(remove)){
    chrRemove <- c(chrRemove, which(seqNames %in% remove))
  }
  if(length(chrRemove) > 0){
    chrKeep <- seqNames[-chrRemove]
  }else{
    chrKeep <- seqNames
  }
  #this function restores seqlevels
  seqlevels(gr, pruning.mode=pruningMode) <- chrKeep
  
  return(gr)
  
}

#' Retreive a non-overlapping set of regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be used to determine how overlapping regions should be resolved.
#' The resolution of overlapping regions also depends on `decreasing`. For example, if a column named "score" is used for `by`,
#' `decreasing = TRUE` means that the highest "score" in the overlap will be retained and `decreasing = FALSE` means that the
#' lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the values in the column indicated via `by` should be ordered in decreasing
#' order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose A boolean value indicating whether the output should include extra reporting.
#' @export
nonOverlappingGR <- function(
    gr = NULL, 
    by = "score", 
    decreasing = TRUE, 
    verbose = FALSE
){
  
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = by, name = "by", valid = c("character"))
  .validInput(input = decreasing, name = "decreasing", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  
  stopifnot(by %in% colnames(mcols(gr)))
  
  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", decreasing = TRUE){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  i <-  0
  grConverge <- gr
  while(length(grConverge) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i <-  i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge, 
      filter = TRUE, 
      by = by, 
      decreasing = decreasing)
    
    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect, 
      invert=TRUE, 
      ignore.strand = TRUE) #blacklist selected gr
    
    if(i == 1){ #if i=1 then set gr_all to clustered
      grAll <- grSelect
      
    }else{
      grAll <- c(grAll, grSelect)
    } 
    
  }
  message(sprintf("Converged after %s iterations!", i))
  
  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))
  
  return(grAll)
  
}

#' Extend regions from a Genomic Ranges object
#'
#' This function extends each region in a Genomic Ranges object by a designated upstream and downstream extension in a strand-aware fashion
#'
#' @param gr A `GRanges` object.
#' @param upstream The number of basepairs upstream (5') to extend each region in `gr` in a strand-aware fashion.
#' @param downstream The number of basepairs downstream (3') to extend each region in `gr` in a strand-aware fashion.
#' @export
extendGR <-  function(gr = NULL, upstream = NULL, downstream = NULL){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
  isMinus <- BiocGenerics::which(strand(gr) == "-")
  isOther <- BiocGenerics::which(strand(gr) != "-")
  #Forward
  start(gr)[isOther] <- start(gr)[isOther] - upstream
  end(gr)[isOther] <- end(gr)[isOther] + downstream
  #Reverse
  end(gr)[isMinus] <- end(gr)[isMinus] + upstream
  start(gr)[isMinus] <- start(gr)[isMinus] - downstream
  return(gr)
}

##########################################################################################
# ggplot2 Wrapper Methods For Easy Plotting
##########################################################################################

#' A ggplot-based dot plot wrapper function
#'
#' This function is a wrapper around ggplot geom_point to allow for a more intuitive plotting of ArchR data.
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param color A numeric/categorical vector used to determine the coloration for each point.
#' @param discrete A boolean value indicating whether the supplied data is discrete (`TRUE`) or continuous (`FALSE`).
#' @param discreteSet The name of a custom palette from `ArchRPalettes` to use for categorical/discrete color.
#' This argument is only used if `discrete` is set to `TRUE`.
#' @param continuousSet The name of a custom palette from `ArchRPalettes` to use for numeric color.
#' This argument is only used if `discrete` is set to `FALSE`.
#' @param labelMeans A boolean value indicating whether the mean of each categorical/discrete color should be labeled.
#' @param pal A custom palette used to override discreteSet/continuousSet for coloring vector.
#' @param defaultColor The default color for points that do not have another color applied (i.e. `NA` values).
#' @param highlightPoints A integer vector describing which points to hightlight. The remainder of points will be colored light gray.
#' @param colorDensity A boolean value indicating whether the density of points on the plot should be indicated by color.
#' If `TRUE`, continuousSet is used as the color palette.
#' @param size The numeric size of the points to be plotted.
#' @param xlim A numeric vector of two values indicating the lower and upper bounds of the x-axis on the plot.
#' @param ylim A numeric vector of two values indicating the lower and upper bounds of the y-axis on the plot.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum and minimum
#' values if `xlim` and `ylim` are not provided. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param randomize A boolean value indicating whether to randomize the order of the points when plotting.
#' @param seed A numeric seed number for use in randomization.
#' @param colorTitle A title to be added to the legend if `color` is supplied.
#' @param colorOrder A vector that allows you to control the order of palette colors associated with the values in `color`.
#' For example if you have `color` as `c("a","b","c")` and want to have the first color selected from the palette be used for
#' "c", the second color for "b", and the third color for "a", you would supply the `colorOrder` as `c("c", "b", "a")`.
#' @param colorLimits A numeric vector of two values indicating the lower and upper bounds of colors if numeric. Values
#' beyond these limits are thresholded.
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param legendSize The size in inches to use for plotting the color legend.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param labelAsFactors A boolean indicating whether to label the `color` input as a numeric factor (`TRUE`) or with a character string (`FALSE`).
#' @param fgColor The foreground color of the plot.
#' @param bgColor The background color of the plot.
#' @param bgWidth The background relative width size of the halos in the labeling.
#' @param labelSize The numeric font size of labels.
#' @param addFit A string indicating the method to use for adding a fit/regression line to the plot (see `ggplot2::geom_smooth()` methods).
#' If set to `NULL`, no fit/regression line is added.
#' @param rastr A boolean value that indicates whether the plot should be rasterized using `ggrastr`. This does not rasterize
#' lines and labels, just the internal portions of the plot.
#' @param dpi The resolution in dots per inch to use for the plot.
#' @export
ggPoint <- function(
    x = NULL, 
    y = NULL, 
    color = NULL, 
    discrete = TRUE, 
    discreteSet = "stallion",
    continuousSet = "solarExtra", 
    labelMeans = TRUE,  
    pal = NULL, 
    defaultColor = "lightGrey",
    highlightPoints = NULL,
    colorDensity = FALSE,
    size = 1, 
    xlim = NULL, 
    ylim = NULL, 
    extend = 0.05, 
    xlabel = "x", 
    ylabel = "y", 
    title = "", 
    randomize = FALSE, 
    seed = 1,
    colorTitle = NULL, 
    colorOrder = NULL, 
    colorLimits = NULL,
    alpha = 1, 
    baseSize = 10, 
    legendSize = 3,
    ratioYX = 1, 
    labelAsFactors = TRUE,
    fgColor = "black", 
    bgColor = "white", 
    bgWidth = 1,
    labelSize = 3,
    addFit = NULL, 
    rastr = FALSE, 
    dpi = 300,
    ...
){
  
  .validInput(input = x, name = "x", valid = c("numeric"))
  .validInput(input = y, name = "y", valid = c("numeric"))
  .validInput(input = color, name = "color", valid = c("numeric", "character", "null"))
  .validInput(input = discrete, name = "discrete", valid = c("boolean"))
  .validInput(input = discreteSet, name = "discreteSet", valid = c("character"))
  .validInput(input = continuousSet, name = "continuousSet", valid = c("character"))
  .validInput(input = labelMeans, name = "labelMeans", valid = c("boolean"))
  .validInput(input = pal, name = "pal", valid = c("character", "null"))
  .validInput(input = defaultColor, name = "defaultColor", valid = c("character"))
  .validInput(input = highlightPoints, name = "highlightPoints", valid = c("integer", "null"))
  .validInput(input = colorDensity, name = "colorDensity", valid = c("boolean"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = xlim, name = "xlim", valid = c("numeric", "null"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = extend, name = "extend", valid = c("numeric"))
  .validInput(input = xlabel, name = "xlabel", valid = c("character"))
  .validInput(input = ylabel, name = "ylabel", valid = c("character"))
  .validInput(input = title, name = "title", valid = c("character"))
  .validInput(input = randomize, name = "randomize", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = colorTitle, name = "colorTitle", valid = c("character", "null"))
  .validInput(input = colorOrder, name = "colorOrder", valid = c("character", "null"))
  .validInput(input = colorLimits, name = "colorLimits", valid = c("numeric", "null"))
  .validInput(input = alpha, name = "alpha", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = legendSize, name = "legendSize", valid = c("numeric"))
  .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric"))
  .validInput(input = labelAsFactors, name = "labelAsFactors", valid = c("boolean"))
  .validInput(input = fgColor, name = "fgColor", valid = c("character", "null"))
  .validInput(input = bgColor, name = "bgColor", valid = c("character"))
  .validInput(input = bgWidth, name = "bgWidth", valid = c("numeric"))
  .validInput(input = labelSize, name = "labelSize", valid = c("numeric"))
  .validInput(input = addFit, name = "addFit", valid = c("character", "null"))
  .validInput(input = rastr, name = "rastr", valid = c("boolean"))
  .validInput(input = dpi, name = "dpi", valid = c("numeric"))
  
  stopifnot(length(y) == length(x))
  if(length(x) < 5){
    stop("x must be at least length 5 to plot!")
  }
  
  if(randomize){
    set.seed(seed)
    idx <- sample(seq_along(x), length(x))
  }else{
    idx <- seq_along(x)
  }
  
  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    x <- x[include]
    y <- y[include]
    if(!is.null(color)){
      color <- color[include]
    }
  }
  
  if(is.null(xlim)){
    xlim <- range(df$x) %>% extendrange(f = extend)
  }
  
  if(is.null(ylim)){
    ylim <- range(df$y) %>% extendrange(f = extend)
  }
  
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  
  #Plot
  .requirePackage("ggplot2", source = "cran")
  
  if (is.null(color) & !colorDensity) {
    
    p <- ggplot(df[idx,], aes(x = x, y = y)) + 
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
      xlab(xlabel) + ylab(ylabel) + 
      ggtitle(title) + 
      theme_ArchR(baseSize = baseSize)
    
    if(rastr){
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
      # if(!requireNamespace("ggrastr", quietly = TRUE)){
      #   message("ggrastr is not available for rastr of points, continuing without rastr!")
      #   p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
      # }else{
      #   .requirePackage("ggrastr")
      #   p <- p + geom_point_rast(
      #       size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
      # }
    }else{
      p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
    }
    
  }else {
    
    if(colorDensity){
      
      discrete <- FALSE
      df <- .getDensity(x, y, n = 100, sample = NULL) #change
      df <- df[order(df$density), ,drop=FALSE]
      df$color <- df$density
      
      if(is.null(colorTitle)){
        colorTitle <- "density"
      }
      
    }else if(discrete){
      
      if(!is.null(highlightPoints)){
        if(length(highlightPoints) < length(color)){
          color[-highlightPoints] <- "Non.Highlighted"
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      color <- paste0(color)
      
      if(!is.null(colorOrder)){
        if(!all(color %in% colorOrder)){
          stop("Not all colors are in colorOrder!")
        }
      }else{
        colorOrder <- gtools::mixedsort(unique(color))
      }
      
      if(is.null(colorTitle)){
        colorTitle <- "color"
      }
      
      stopifnot(length(color) == nrow(df))
      df$color <- factor(color, levels = colorOrder)
      
      if(labelAsFactors){
        df$color <- factor(
          x = paste0(paste0(match(paste0(df$color), paste0(levels(df$color)))), "-", paste0(df$color)), 
          levels = paste0(seq_along(levels(df$color)), "-", levels(df$color))
        )
        if(!is.null(pal)){
          #print(pal)
          #print(paste0(levels(df$color))[match(names(pal), colorOrder)])
          names(pal) <- paste0(levels(df$color))[match(names(pal), colorOrder)]
        }
        colorOrder <- paste0(levels(df$color))
      }
      
    }else{
      stopifnot(length(color) == nrow(df))
      if(!is.null(highlightPoints)){
        if(length(highlightPoints) < length(color)){
          color[-highlightPoints] <- NA
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      if(!is.null(colorLimits)){
        color[color < min(colorLimits)] <- min(colorLimits)
        color[color > max(colorLimits)] <- max(colorLimits)
      }
      df$color <- color
    }
    
    p <- ggplot(df[idx,], aes(x = x, y = y, color = color)) +  
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) + 
      xlab(xlabel) + ylab(ylabel) + 
      ggtitle(title) + theme_ArchR(baseSize = baseSize) +
      theme(legend.direction = "horizontal", legend.box.background = element_rect(color = NA)) +
      labs(color = colorTitle)
    
    if(rastr){
      
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha, 
        raster.width = min(par('fin')), 
        raster.height = (ratioYX * min(par('fin')))
      )
      
      # if(!requireNamespace("ggrastr", quietly = TRUE)){
      #   message("ggrastr is not available for rastr of points, continuing without rastr!")
      #   message("To install ggrastr try : devtools::install_github('VPetukhov/ggrastr')")
      #   p <- p + geom_point(size = size, alpha = alpha)
      # }else{
      #   .requirePackage("ggrastr", installInfo = "devtools::install_github('VPetukhov/ggrastr')")
      #   p <- p + geom_point_rast(
      #       size = size, raster.dpi = dpi, alpha = alpha, 
      #       raster.width=par('fin')[1], 
      #       raster.height = (ratioYX * par('fin')[2])
      #     )
      # }
      
    }else{
      
      p <- p + geom_point(size = size, alpha = alpha)
      
    }
    
    if (discrete) {
      
      if (!is.null(pal)) {
        p <- p + scale_color_manual(values = pal)
      }else {
        pal <- paletteDiscrete(set = discreteSet, values = colorOrder)
        if(!is.null(highlightPoints)){
          pal[grep("Non.Highlighted", names(pal))] <- "lightgrey"
        }
        #print(pal)
        p <- p + scale_color_manual(values = pal) +
          guides(color = guide_legend(override.aes = list(size = legendSize, shape = 15)))
      }
      
      if (labelMeans) {
        
        dfMean <- split(df, df$color) %>% lapply(., function(x) {
          data.frame(x = median(x[, 1]), y = median(x[, 2]), color = x[1, 3])
        }) %>% Reduce("rbind", .)
        
        if(labelAsFactors){
          dfMean$label <- stringr::str_split(paste0(seq_len(nrow(dfMean))), pattern = "\\-", simplify=TRUE)[,1]
        }else{
          dfMean$label <- dfMean$color
        }
        dfMean$text <- stringr::str_split(dfMean$color, pattern = "-", simplify = TRUE)[,1]
        
        # make halo layers, similar to https://github.com/GuangchuangYu/shadowtext/blob/master/R/shadowtext-grob.R#L43
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        xo <- bgWidth * diff(range(df$x)) / 300
        yo <- bgWidth * diff(range(df$y)) / 300
        for (i in theta) {
          p <- p + 
            geom_text(data = dfMean, 
                      aes_q(
                        x = bquote(x + .(cos(i) * xo)),
                        y = bquote(y + .(sin(i) * yo)),
                        label = ~text
                      ),
                      size = labelSize,
                      color = bgColor
            )
        }
        
        if(is.null(fgColor)){
          p <- p + geom_text(data = dfMean, aes(x = x, y = y, color = color, label = label), size = labelSize, show.legend = FALSE)
        }else{
          p <- p + geom_text(data = dfMean, aes(x = x, y = y, label = label), color = fgColor, size = labelSize, show.legend = FALSE) 
        }
        
      }
      
    }else{
      
      if (!is.null(pal)) {
        if(!is.null(colorLimits)){
          p <- p + scale_colour_gradientn(colors = pal, limits=colorLimits, na.value = "lightgrey")
        }else{
          p <- p + scale_colour_gradientn(colors = pal, na.value = "lightgrey")
        }
      }else {
        if(!is.null(colorLimits)){
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), limits=colorLimits, na.value = "lightgrey")
        }else{
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), na.value = "lightgrey")
        }
      }
    }
    
  }
  
  if (!is.null(addFit)) {
    p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, color = "black") + 
      ggtitle(paste0(title, "\nPearson = ", round(cor(df$x, df$y), 3), "\nSpearman = ", round(cor(df$x, df$y, method = "spearman"), 3)))
  }
  
  p <- p + theme(legend.position = "bottom", legend.key = element_rect(size = 2))#, legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(0.1, 'cm'))
  
  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }
  
  return(p)
  
}

#' A ggplot-based one-to-one dot plot wrapper function
#'
#' This function is a wrapper around ggplot geom_point to allow for plotting one-to-one sample comparisons in ArchR.
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param size The numeric size of the points to plot.
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param min The lower limit of the x and y axes as a numeric quantile between 0 and 1.
#' @param max The upper limit of the x and y axes as a numeric quantile between 0 and 1.
#' @param nPlot The number of points to plot. When this value is less than the total points, the `sample` function is used to extract random data points to be plotted.
#' @param nKernel The number of grid points in each direction to use when computing the kernel with `MASS::kde2d()`.
#' @param densityMax The quantile that should be represented by the maximum color on the continuous scale designated by `pal`. Values above `densityMax` will be thresholded to the maximum color on the color scale.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum value on either axis. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end beyond `quantile(c(x,y), max)` and `quantile(c(x,y), min)`.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels, just the internal portions of the plot.
#' @param pal A custom palette from `ArchRPalettes` used to display the density of points on the plot.
#' @param ... Additional params to be supplied to ggPoint
#' @export
ggOneToOne <- function (
    x = NULL,
    y = NULL,
    size = 2, 
    alpha = 1,
    xlabel = "x", 
    ylabel = "y", 
    title = "Correlation",
    min = 0.05, 
    max = 0.9999, 
    nPlot = 100 * 10^3, 
    nKernel = 100, 
    densityMax = 0.95, 
    extend = 0.05, 
    baseSize = 6, 
    rastr = TRUE,
    pal = paletteContinuous(set = "blueYellow"),
    ...
){
  
  .validInput(input = x, name = "x", valid = c("numeric"))
  .validInput(input = y, name = "y", valid = c("numeric"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = alpha, name = "alpha", valid = c("numeric"))
  .validInput(input = xlabel, name = "xlabel", valid = c("character"))
  .validInput(input = ylabel, name = "ylabel", valid = c("character"))
  .validInput(input = title, name = "title", valid = c("character"))
  .validInput(input = min, name = "min", valid = c("numeric"))
  .validInput(input = max, name = "max", valid = c("numeric"))
  .validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  .validInput(input = nKernel, name = "nKernel", valid = c("numeric"))
  .validInput(input = densityMax, name = "densityMax", valid = c("numeric"))
  .validInput(input = extend, name = "extend", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = rastr, name = "rastr", valid = c("boolean"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  
  #Check for NA
  idx <- which(!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y))
  x <- x[idx]
  y <- y[idx]
  
  #Ratio X/Y
  lim <- quantile(c(x, y), c(min, max)) %>% extendrange(f = extend)
  ratioXY <- diff(lim)/diff(lim)
  
  #Calculate Correlations
  pearson <- round(cor(x, y, method = "pearson", use = "complete"), 3)
  spearman <- round(cor(x, y, method = "spearman", use = "complete"), 3)
  title <- sprintf("%s \nPearson = %s , Spearman = %s", title, pearson, spearman)
  
  #Get Density
  message("adding denisty..")
  df <- .getDensity(x, y, n = nKernel, sample = nPlot) #change
  df <- df[order(df[, "density"]), ]
  
  #GGPlot
  message("plotting..")
  gg <- ggPoint(
    x = df$x, 
    y = df$y, 
    color = df$density, 
    pal = pal,
    xlabel = xlabel,
    ylabel = ylabel,
    discrete = FALSE, 
    colorTitle = "density",
    xlim = lim, 
    ylim = lim, 
    size = size, 
    alpha = alpha, 
    title = title, 
    baseSize = baseSize,
    rastr = rastr,
    ...
  ) + geom_abline(slope = 1, intercept = 0, lty = "dashed")
  
  return(gg)
  
}

.getDensity <- function(x = NULL, y = NULL, n = 100, sample = NULL, densityMax = 0.95){
  #modified from http://slowkow.com/notes/ggplot2-color-by-density/
  df <- data.frame(x=x,y=y)
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  df$density <- dens$z[ii]
  df$density[df$density > quantile(unique(df$density),densityMax)] <- quantile(unique(df$density),densityMax) #make sure the higher end doesnt bias colors
  if(!is.null(sample)){
    df <- df[sample(nrow(df), min(sample,nrow(df))),]
  }
  return(df)
}

#' A ggplot-based Hexplot wrapper function summary of points in a standardized manner
#'
#' This function will plot x,y coordinate values summarized in hexagons in a standardized manner
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param color A numeric/categorical vector containing coloring information for each point.
#' @param pal A custom continuous palette from `ArchRPalettes` for coloration of hexes.
#' @param bins The number of bins to be used for plotting the hexplot. `bins` indicates the total number of hexagons that will fit within the surface area of the plot. 
#' @param xlim A numeric vector of two values indicating the lower and upper bounds of the x-axis on the plot.
#' @param ylim A numeric vector of two values indicating the lower and upper bounds of the y-axis on the plot.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum and minimum values if `xlim` and `ylim` are not provided. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param colorTitle The label to use for the legend corresponding to `color`.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param FUN The function to use for summarizing data into hexagons. Typically "mean" or something similar.
#' @param hexCut If this is not null, a quantile cut is performed to threshold the top and bottom of the distribution of values.
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(a,b) where `a` is the upper threshold
#' and `b` is the lower threshold. For example, hexCut = c(0.025,0.975) will take the top and bottom 2.5 percent of values and set
#' them to the value of the 97.5th and 2.5th percentile values respectively.
#' @param addPoints A boolean value indicating whether individual points should be shown on the hexplot.
#' @param ... Additional params for plotting
#' @export
ggHex <- function(
    x = NULL, 
    y = NULL, 
    color = NULL, 
    pal = paletteContinuous(set = "solarExtra"), 
    bins = 200,
    xlim = NULL, 
    ylim = NULL, 
    extend = 0.05, 
    xlabel = "x", 
    ylabel = "y",
    title = "", 
    colorTitle = "values", 
    baseSize = 6,
    ratioYX = 1, 
    FUN = "median", 
    hexCut = c(0.02, 0.98),
    addPoints = FALSE,
    ...
){
  
  .validInput(input = x, name = "x", valid = c("numeric"))
  .validInput(input = y, name = "y", valid = c("numeric"))
  .validInput(input = color, name = "color", valid = c("numeric"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  .validInput(input = bins, name = "bins", valid = c("integer"))
  .validInput(input = xlim, name = "xlim", valid = c("numeric", "null"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = xlabel, name = "xlabel", valid = c("character"))
  .validInput(input = ylabel, name = "ylabel", valid = c("character"))
  .validInput(input = title, name = "title", valid = c("character"))
  .validInput(input = colorTitle, name = "colorTitle", valid = c("character", "null"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric"))
  .validInput(input = FUN, name = "FUN", valid = c("character"))
  .validInput(input = hexCut, name = "quantCut", valid = c("numeric", "null"))
  .validInput(input = addPoints, name = "addPoints", valid = c("boolean"))
  
  #require hexbin to be installed. otherwise, this section wont work properly
  .requirePackage(x = "hexbin", source = "CRAN")
  
  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    if(!is.null(color)){
      color <- color[include]
    }
  }
  df$color <- color
  
  if (is.null(xlim)) {
    xlim <- range(df$x) %>% extendrange(f = extend)
  }
  if (is.null(ylim)) {
    ylim <- range(df$y) %>% extendrange(f = extend)
  }
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  
  p <- ggplot()
  
  if(addPoints){
    p <- p + .geom_point_rast2(data = df, aes(x=x,y=y), color = "lightgrey")
    # if(requireNamespace("ggrastr", quietly = TRUE)){
    #   .requirePackage("ggrastr", installInfo = "devtools::install_github('VPetukhov/ggrastr')")
    #   p <- p + geom_point_rast(data = df, aes(x=x,y=y), color = "lightgrey")
    # }else{
    #   message("ggrastr is not available for rastr of points, continuing without points!")
    #   message("To install ggrastr try : devtools::install_github('VPetukhov/ggrastr')")
    #}
  }
  
  values <- ggplot_build(p + stat_summary_hex(data = df, aes(x=x,y=y,z=color), fun = FUN, bins = bins, color = NA))$data[[1]]$value
  if(!is.null(hexCut)){
    limits <- quantile(values, c(min(hexCut), max(hexCut)), na.rm=TRUE)
  }else{
    limits <- c(min(values), max(values))
  }
  
  p <- p + stat_summary_hex(data = df, aes(x=x,y=y,z=color), fun = FUN, bins = bins, color = NA) +
    scale_fill_gradientn(
      colors = pal,
      limits = limits, 
      oob = scales::squish
    ) +
    xlab(xlabel) + 
    ylab(ylabel) + 
    ggtitle(title) +
    theme_ArchR(baseSize = baseSize) +
    coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) +
    theme(legend.direction="horizontal", legend.box.background = element_rect(color = NA)) +
    labs(fill = colorTitle)
  
  p <- p + theme(legend.position = "bottom")
  
  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }
  
  p
  
}

#' A ggplot-based ridge/violin plot wrapper function
#'
#' This function is a wrapper around ggplot geom_density_ridges or geom_violin to allow for plotting group distribution plots in ArchR.
#' 
#' @param x A character vector containing the categorical x-axis values for each y-axis value.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param groupOrder A character vector indicating a custom order for plotting x-axis categorical values. Should contain all possible
#' values of `x` in the desired order.
#' @param groupSort A boolean indicating whether to sort groups based on the average value of the group.
#' @param size The line width for boxplot/summary lines.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param ridgeScale A numeric indicating the relative size for each ridge in the ridgeplot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param title The title of the plot.
#' @param pal A named custom palette (see `paletteDiscrete()` and `ArchRPalettes`) for discrete coloring.
#' @param addBoxPlot A boolean indicating whether to add a boxplot to the plot if `plotAs="violin"`.
#' @param plotAs A string indicating how the groups should be plotted. Acceptable values are "ridges" (for a `ggrides`-style plot) or "violin" (for a violin plot).
#' @param ... Additional parameters to pass to `ggplot2` for plotting.
#' @export
ggGroup <- function(
    x = NULL, 
    y = NULL, 
    xlabel = NULL, 
    ylabel = NULL, 
    groupOrder = NULL,
    groupSort = FALSE,
    size = 1,  
    baseSize = 10,
    ridgeScale = 1, 
    ratioYX = NULL,
    alpha = 1,
    title = "", 
    pal = paletteDiscrete(values=x, set = "stallion"),
    addBoxPlot = TRUE,
    plotAs = "ridges",
    ...
){
  
  .validInput(input = x, name = "x", valid = c("character"))
  .validInput(input = y, name = "y", valid = c("numeric"))
  .validInput(input = xlabel, name = "xlabel", valid = c("character", "null"))
  .validInput(input = ylabel, name = "ylabel", valid = c("character", "null"))
  .validInput(input = groupOrder, name = "groupOrder", valid = c("character", "null"))
  .validInput(input = groupSort, name = "groupSort", valid = c("boolean"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = ridgeScale, name = "ridgeScale", valid = c("numeric"))
  .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric", "null"))
  .validInput(input = alpha, name = "alpha", valid = c("numeric"))
  .validInput(input = title, name = "title", valid = c("character"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  .validInput(input = addBoxPlot, name = "addBoxPlot", valid = c("boolean"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character"))
  
  names(y) <- x
  dm <- stats::aggregate(y ~ names(y), FUN = mean)
  df <- data.frame(x, y)
  
  if(!is.null(groupOrder)){
    if(!all(x %in% groupOrder)){
      stop("Not all x values are present in groupOrder!")
    }
  }else{
    if(groupSort){
      groupOrder <- paste0(dm[,1])[order(dm[,2], decreasing= FALSE)]
    }else{
      if(tolower(plotAs) == "ridges"){
        groupOrder <- rev(gtools::mixedsort(unique(x)))
      }else{
        groupOrder <- gtools::mixedsort(unique(x))
      }
    }
  }
  
  df$x <- factor(df$x, groupOrder)
  
  p <- ggplot(df, aes(x = x, y = y, color = x)) +
    scale_color_manual(values = pal, guide = "none") + 
    scale_fill_manual(values = pal, guide = "none") +
    ggtitle(title)
  
  if(tolower(plotAs) == "ridges" | tolower(plotAs) == "ggridges"){
    if(!requireNamespace("ggridges", quietly = TRUE)){
      type <- "violin"
      message("ggridges is not available for plotting, continuing with geom_violin!")
      message("To install ggridges try : install.packages('ggridges')")
      p <- p + geom_violin(aes_string(fill="x"), alpha = alpha)
    }else{
      type <- "ridges"
      .requirePackage("ggridges", source = "cran")
      #p <- p + 
      #  stat_density_ridges(aes_string(x = "y", y = "x", fill = "x"), 
      #    quantile_lines = TRUE, quantiles = c(0.5), alpha = alpha, color = "black",
      #    scale = ridgeScale
      #  ) + scale_y_discrete(expand = c(0, 0))
      #   stat_density_ridges(
      #     aes_string(x = "y", y = "x", fill = "x"),
      #     quantile_lines = TRUE,
      #     alpha = alpha,
      #     geom = "density_ridges_gradient",
      #     calc_ecdf = TRUE,
      #     quantiles = c(0.5)
      # )
      val <- 1/length(unique(x))
      p <- p + geom_density_ridges(data = df,
                                   aes(x = y, y = x, color = x, fill = x), scale = ridgeScale,
                                   alpha = alpha, color = "black") + scale_y_discrete(expand = expansion(mult = c(0.01, val)))
    }
  }else{
    type <- "violin"
    p <- p + geom_violin(aes_string(x = "x", y = "y", color = "x", fill="x"), alpha = alpha)
  }
  
  if(addBoxPlot & type == "violin"){
    p <- p + geom_boxplot(size = size, outlier.size = 0, outlier.stroke = 0, fill = NA) 
  }
  
  if(type != "violin"){
    p <- p + theme_ArchR(baseSize = baseSize)
  }else{
    p <- p + theme_ArchR(xText90 = TRUE, baseSize = baseSize)
  }
  
  if(!is.null(ratioYX)){
    p <- p + coord_fixed(ratioYX, expand = TRUE)
  }
  
  if (!is.null(xlabel)) {
    if(type=="violin"){
      p <- p + xlab(xlabel)
    }else{
      p <- p + xlab(ylabel)
    }
  }
  
  if (!is.null(ylabel)) {
    if(type=="violin"){
      p <- p + ylab(ylabel)
    }else{
      p <- p + ylab(xlabel)
    }
  }
  
  p <- p + theme(legend.position = "bottom")
  
  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }
  
  return(p)
  
}

#' Align ggplot plots vertically or horizontally
#'
#' This function aligns ggplots vertically or horizontally
#'
#' @param ... All additional arguments will be interpreted as `ggplot2` plot objects and used if and only if `plotList` is `NULL`
#' @param plotList A list of `ggplot2` plot objects to be aligned.
#' @param sizes A numeric vector or list of values indicating the relative size for each of the objects in `plotList` or supplied in `...`. If the plot is supplied in `...` the order is the same as the input in this function. If set to NULL all plots will be evenly distributed.
#' @param type A string indicating wheter vertical ("v") or horizontal ("h") alignment should be used for the multi-plot layout.
#' @param draw A boolean value indicating whether to draw the plot(s) (`TRUE`) or return a graphical object (`FALSE`).
#' @export
ggAlignPlots <- function(
    ..., 
    plotList = NULL, 
    sizes = NULL, 
    type = "v",  
    draw = TRUE
){
  
  .validInput(input = plotList, name = "plotList", valid = c("list", "null"))
  .validInput(input = sizes, name = "sizes", valid = c("numeric", "null"))
  .validInput(input = type, name = "type", valid = c("character"))
  .validInput(input = draw, name = "draw", valid = c("boolean"))
  if(type %ni% c("v", "h")){
    stop("type must be v (vertical) or h (horizontal)!")
  }
  
  #http://stackoverflow.com/a/21503904
  
  .requirePackage("gtable", source = "cran")
  
  if(is.null(plotList)){
    plotList <- list(...)
  }
  
  ## test that only passing plots
  stopifnot(do.call(all, lapply(plotList, inherits, "gg")))
  
  gl <- lapply(plotList, ggplotGrob)
  
  #if ncols do not match fill with empty gtables_add_cols
  if(type == "v" | type == "vertical"){
    maxCol <- max(unlist(lapply(gl, ncol)))
    gl <- lapply(gl, function(x){
      while(ncol(x) < max(maxCol)){
        x <- gtable::gtable_add_cols(x, unit(1, "null"))
      }
      return(x)
    })
  }
  
  combined <- Reduce(function(x, y)
    if(type == "v" | type == "vertical"){
      gtable:::rbind_gtable(x,y,"first")
    }else{
      gtable:::cbind_gtable(x,y,"first")
    }, gl[-1], gl[[1]])
  
  if(type == "v" | type == "vertical"){
    combined$widths <- do.call(grid::unit.pmax, lapply(gl, "[[", "widths"))
    #remove vertical spaces from background layout
    combined$heights[combined$layout$t[grepl("background", combined$layout$name)][-1]] <- grid::unit(rep(0,length(combined$heights[combined$layout$t[grepl("background", combined$layout$name)][-1]])), "cm")
    if(!missing(sizes)){
      sList <- lapply(seq_along(gl), function(x){
        orig <- gl[[x]]$heights[gl[[x]]$layout$t[grepl("panel", gl[[x]]$layout$name)]]
        new <- rep(sizes[[x]]/length(orig),length(orig))
        return(new)
      })
      s <- grid::unit(unlist(sList), "null")
      combined$heights[combined$layout$t[grepl("panel", combined$layout$name)]] <- s
    }
  }else if(type == "h" | type == "horizontal"){
    combined$heights <- do.call(grid::unit.pmax, lapply(gl, "[[", "heights"))
    if(!missing(sizes)){
      sList <- lapply(seq_along(gl), function(x){
        orig <- gl[[x]]$widths[gl[[x]]$layout$l[grepl("panel", gl[[x]]$layout$name)]]
        new <- rep(sizes[[x]]/length(orig),length(orig))
        return(new)
      })
      s <- grid::unit(unlist(sList), "null")
      combined$widths[combined$layout$l[grepl("panel", combined$layout$name)]] <- s
    }
  }else{
    stop("Unrecognized type ", type)
  }
  
  if(draw){
    grid::grid.newpage()
    grid::grid.draw(combined)
  }else{
    combined    
  }
  
}

#' ggplot2 default theme for ArchR
#'
#' This function returns a ggplot2 theme that is black borded with black font.
#' 
#' @param color The color to be used for text, lines, ticks, etc for the plot.
#' @param textFamily The font default family to be used for the plot.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param baseLineSize The base line width (in points) to be used throughout the plot.
#' @param baseRectSize The base line width (in points) to use for rectangular boxes throughout the plot.
#' @param plotMarginCm The width in centimeters of the whitespace margin around the plot.
#' @param legendPosition The location to put the legend. Valid options are "bottom", "top", "left", and "right.
#' @param legendTextSize The base text size (in points) for the legend text.
#' @param axisTickCm The length in centimeters to be used for the axis ticks.
#' @param xText90 A boolean value indicating whether the x-axis text should be rotated 90 degrees counterclockwise.
#' @param yText90 A boolean value indicating whether the y-axis text should be rotated 90 degrees counterclockwise.
#' @export
theme_ArchR <- function(
    color = "black",
    textFamily = "sans",
    baseSize = 10, 
    baseLineSize = 0.5,
    baseRectSize = 0.5,
    plotMarginCm = 1,
    legendPosition = "bottom",
    legendTextSize = 5,
    axisTickCm = 0.1,
    xText90 = FALSE,
    yText90 = FALSE
){
  
  .validInput(input = color, name = "color", valid = c("character"))
  .validInput(input = textFamily, name = "textFamily", valid = c("character"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = baseLineSize, name = "baseLineSize", valid = c("numeric"))
  .validInput(input = baseRectSize, name = "baseRectSize", valid = c("numeric"))
  .validInput(input = plotMarginCm, name = "plotMarginCm", valid = c("numeric"))
  .validInput(input = legendPosition, name = "legendPosition", valid = c("character"))
  .validInput(input = legendTextSize, name = "legendTextSize", valid = c("numeric"))
  .validInput(input = axisTickCm, name = "axisTickCm", valid = c("numeric"))
  .validInput(input = xText90, name = "xText90", valid = c("boolean"))
  .validInput(input = yText90, name = "yText90", valid = c("boolean"))
  
  theme <- theme_bw() + theme(
    text = element_text(family = textFamily),
    axis.text = element_text(color = color, size = baseSize), 
    axis.title = element_text(color = color, size = baseSize),
    title = element_text(color = color, size = baseSize),
    plot.margin = unit(c(plotMarginCm, plotMarginCm, plotMarginCm, plotMarginCm), "cm"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = color, size = (4/3) * baseRectSize * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    axis.ticks.length = unit(axisTickCm, "cm"), 
    axis.ticks = element_line(color = color, size = baseLineSize * (4/3) * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(color = color, size = legendTextSize),
    legend.box.background = element_rect(color = NA),
    #legend.box.background = element_rect(fill = "transparent"),
    legend.position = legendPosition,
    strip.text = element_text(size = baseSize, color="black")#,
    #plot.background = element_rect(fill = "transparent", color = NA)
  )
  
  if(xText90){
    theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  if(yText90){
    theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90, vjust = 1))
  }
  
  return(theme)
  
}



##########################################################################################
# ggplot2 helper functions
##########################################################################################

.checkCairo <- function(){
  tryCatch({
    tmp <- dev.cur()
    Cairo::Cairo(type='raster')
    dev.off()
    dev.set(tmp)
    TRUE
  }, error = function(e){
    FALSE
  })
}

## Adapted from 
## https://github.com/tidyverse/ggplot2/blob/660aad2db2b3495ae0d8040915a40d247133ffc0/R/geom-point.r
## from https://github.com/VPetukhov/ggrastr/blob/master/R/geom-point-rast.R
## This funciton now handles issues with Cairo installation that can lead to plot errors
.geom_point_rast2 <- function(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    ...,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE,
    raster.width = min(par('fin')), 
    raster.height = min(par('fin')), 
    raster.dpi = 300
){
  
  GeomPointRast <- tryCatch({
    
    if(!.checkCairo()){
      stop()
    }
    
    #Try to create a geom rast for points if not then just use normal geom_point
    ggplot2::ggproto(
      "GeomPointRast",
      ggplot2::GeomPoint,
      required_aes = c("x", "y"),
      non_missing_aes = c("size", "shape", "colour"),
      default_aes = aes(
        shape = 19, colour = "black", size = 1.5, fill = NA,
        alpha = NA, stroke = 0.5
      ),
      
      draw_panel = function(data, panel_params, coord, na.rm = FALSE, 
                            raster.width=min(par('fin')), raster.height=min(par('fin')), raster.dpi=300){
        
        #From ggrastr  
        prevDevID <- dev.cur()
        
        p <- ggplot2::GeomPoint$draw_panel(data, panel_params, coord)
        
        devID <- Cairo::Cairo(
          type='raster', 
          width=raster.width*raster.dpi, 
          height=raster.height*raster.dpi, 
          dpi=raster.dpi, 
          units='px', 
          bg="transparent"
        )[1]
        
        grid::pushViewport(grid::viewport(width=1, height=1))
        
        grid::grid.points(
          x=p$x, 
          y=p$y, 
          pch = p$pch, 
          size = p$size,
          name = p$name, 
          gp = p$gp, 
          vp = p$vp, 
          draw = TRUE
        )
        
        grid::popViewport()
        gridCapture <- grid::grid.cap()
        
        dev.off(devID)
        
        dev.set(prevDevID)
        
        grid::rasterGrob(
          gridCapture, 
          x=0, 
          y=0, 
          width = 1,
          height = 1,
          default.units = "native",
          just = c("left","bottom")
        )
        
      }
      
    )
    
  }, error = function(e){
    
    if(.checkCairo()){
      message("WARNING: Error found with trying to rasterize geom. Continuing without rasterization.")
    }else{
      message("WARNING: Error found with Cairo installation. Continuing without rasterization.")
    }
    
    #Default geom_point
    ggplot2::ggproto(
      "GeomPoint", 
      ggplot2::GeomPoint,
      required_aes = c("x", "y"),
      non_missing_aes = c("size", "shape", "colour"),
      default_aes = aes(
        shape = 19, colour = "black", size = 1.5, fill = NA,
        alpha = NA, stroke = 0.5
      ),
      
      draw_panel = function(data, panel_params, coord, na.rm = FALSE, 
                            raster.width=min(par('fin')), raster.height=min(par('fin')), raster.dpi=300){
        if (is.character(data$shape)) {
          data$shape <- ggplot2:::translate_shape_string(data$shape) #Hidden ggplot2
        }
        
        coords <- coord$transform(data, panel_params)
        
        pGrob <- grid::pointsGrob(
          x = coords$x, 
          y = coords$y,
          pch = coords$shape,
          gp = grid::gpar(
            col = scales::alpha(coords$colour, coords$alpha),
            fill = scales::alpha(coords$fill, coords$alpha),
            # Stroke is added around the outside of the point
            fontsize = coords$size * .pt + coords$stroke * .stroke / 2,
            lwd = coords$stroke * .stroke / 2
          )
        )
        
        pGrob
        
      },
      
      draw_key = ggplot2::draw_key_point
    )
    
    
  })
  
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomPointRast,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      raster.width=raster.width,
      raster.height=raster.height,
      raster.dpi=raster.dpi,
      ...
    )
  )
  
}

####################################################################
# Reading fragments from Arrow Files
####################################################################

#' Get the fragments from an ArchRProject 
#' 
#' This function retrieves the fragments from a given ArchRProject as a GRangesList object.
#'
#' @param ArchRProject An `ArchRProject` object to get fragments from.
#' @param subsetBy A Genomic Ranges object to subset fragments by.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells
#' from the provided ArrowFile using `getCellNames()`.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to `FALSE` for a cleaner output.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getFragmentsFromProject <- function(
    ArchRProj = NULL,
    subsetBy = NULL,
    cellNames = NULL,
    verbose = FALSE,
    logFile = createLogFile("getFragmentsFromProject")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = subsetBy, name = "subsetBy", valid = c("GRanges", "null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  
  if(!is.null(subsetBy)){
    chr <- paste0(unique(seqnames(subsetBy)))
  }else{
    chr <- NULL
  }
  
  ArchR:::.startLogging(logFile = logFile)
  
  FragmentsList <- lapply(seq_along(ArrowFiles), function(x){
    message(sprintf("Reading ArrowFile %s of %s", x, length(ArrowFiles)))
    fragx <- getFragmentsFromArrow(
      ArrowFile = ArrowFiles[x], 
      chr = chr, 
      cellNames = cellNames, 
      verbose = verbose,
      logFile = logFile
    )
    if(!is.null(subsetBy)){
      fragx <- subsetByOverlaps(fragx, subsetBy, ignore.strand = TRUE)
    }
    fragx
  }) %>% SimpleList
  
  names(FragmentsList) <- names(ArrowFiles)
  
  FragmentsList
  
}

#' Get the fragments from an ArrowFile 
#' 
#' This function retrieves the fragments from a given ArrowFile as a GRanges object.
#'
#' @param ArrowFile The path to the ArrowFile from which fragments should be obtained.
#' @param chr A name of a chromosome to be used to subset the fragments `GRanges` object to a specific chromsome if desired.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells
#' from the provided ArrowFile using `getCellNames()`.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to `FALSE` for a cleaner output.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getFragmentsFromArrow <- function(
    ArrowFile = NULL, 
    chr = NULL, 
    cellNames = NULL, 
    verbose = TRUE,
    logFile = createLogFile("getFragmentsFromArrow")
){
  
  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = chr, name = "chr", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getFragmentsFromArrow Input-Parameters", logFile = logFile)
  
  ArrowFile <- .validArrow(ArrowFile)
  
  if(is.null(chr)){
    chr <- .availableSeqnames(ArrowFile, subGroup = "Fragments")
  }
  
  if(any(chr %ni% .availableSeqnames(ArrowFile, subGroup = "Fragments"))){
    stop("Error Chromosome not in ArrowFile!")
  }
  
  out <- lapply(seq_along(chr), function(x){
    .logDiffTime(sprintf("Reading Chr %s of %s", x, length(chr)), t1 = tstart, verbose = verbose, logFile = logFile)
    .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = chr[x], 
      out = "GRanges", 
      cellNames = cellNames, 
      method = "fast"
    )
  })
  
  
  .logDiffTime("Merging", tstart, t1 = tstart, verbose = verbose, logFile = logFile)
  
  out <- tryCatch({
    
    o <- .suppressAll(unlist(GRangesList(out, compress = FALSE)))
    
    if(.isGRList(o)){
      stop("Still a GRangesList")
    }
    
    o
    
  }, error = function(x){
    
    o <- c()
    
    for(i in seq_along(out)){
      if(!is.null(out[[i]])){
        if(i == 1){
          o <- out[[i]]
        }else{
          o <- c(o, out[[i]])
        }
      }
    }
    
    o
    
  })
  
  out
  
}

.getFragsFromArrow <- function(
    ArrowFile = NULL, 
    chr = NULL, 
    out = "GRanges", 
    cellNames = NULL, 
    method = "fast"
){
  
  if(is.null(chr)){
    stop("Need to provide chromosome to read!")
  }
  
  o <- h5closeAll()
  ArrowFile <- .validArrow(ArrowFile)
  
  avSeq <- .availableSeqnames(ArrowFile)
  if(chr %ni% avSeq){
    stop(paste0("Chromosome ", chr ," not in ArrowFile! Available Chromosomes are : ", paste0(avSeq, collapse=",")))
  }
  
  #Get Sample Name
  sampleName <- .h5read(ArrowFile, paste0("Metadata/Sample"), method = method)
  
  o <- h5closeAll()
  nFrags <- sum(.h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method))
  
  if(nFrags==0){
    if(tolower(out)=="granges"){
      output <- GRanges(seqnames = chr, IRanges(start = 1, end = 1), RG = "tmp")
      output <- output[-1,]
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- output[-1,]
    }
    return(output)
  }
  
  if(is.null(cellNames) | tolower(method) == "fast"){
    
    output <- .h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), method = method) %>% 
      {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$RG <- Rle(
      values = paste0(sampleName, "#", .h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues"), method = method)), 
      lengths = .h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method)
    )
    if(!is.null(cellNames)){
      output <- output[BiocGenerics::which(mcols(output)$RG %bcin% cellNames)]
    }
    
  }else{
    
    if(!any(cellNames %in% .availableCells(ArrowFile))){
      
      stop("None of input cellNames are in ArrowFile availableCells!")
      
    }else{
      
      barRle <- Rle(h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues")), h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths")))
      barRle@values <- paste0(sampleName, "#", barRle@values)
      idx <- BiocGenerics::which(barRle %bcin% cellNames)
      if(length(idx) > 0){
        output <- h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), index = list(idx, 1:2)) %>% 
          {IRanges(start = .[,1], width = .[,2])}
        mcols(output)$RG <- barRle[idx]
      }else{
        output <- IRanges(start = 1, end = 1)
        mcols(output)$RG <- c("tmp")
        output <- output[-1,]
      }
    }
    
  }
  
  o <- h5closeAll()
  
  if(tolower(out)=="granges"){
    if(length(output) > 0){
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)    
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)
      output <- output[-1,]
    }
  }
  
  return(output)
}

####################################################################
# Reading Matrices/Arrays from Arrow Files
####################################################################

#' Get a data matrix stored in an ArchRProject
#' 
#' This function gets a given data matrix from an `ArchRProject` and returns it as a `SummarizedExperiment`.
#' This function will return the matrix you ask it for, without altering that matrix unless you tell it to.
#' For example, if you added your `PeakMatrix` using `addPeakMatrix()` with `binarize = TRUE`, then
#' `getMatrixFromProject()` will return a binarized `PeakMatrix`. Alternatively, you could set `binarize = TRUE`
#' in the parameters passed to `getMatrixFromProject()` and the `PeakMatrix` will be binarized as you pull
#' it out. No other normalization is applied to the matrix by this function.
#'
#' @param ArchRProj An `ArchRProject` object to get data matrix from.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return.
#' This is often desired when working with insertion counts. Note that if the matrix has already been binarized previously, this should be set to `TRUE`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromProject <- function(
    ArchRProj = NULL,
    useMatrix = "GeneScoreMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromProject Input-Parameters", logFile = logFile)
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  
  cellNames <- ArchRProj$cellNames
  
  avMat <- getAvailableMatrices(ArchRProj)
  if(useMatrix %ni% avMat){
    stop("useMatrix is not in Available Matrices see getAvailableMatrices")
  }
  
  seL <- .safelapply(seq_along(ArrowFiles), function(x){
    
    .logDiffTime(paste0("Reading ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
                 t1 = tstart, verbose = FALSE, logFile = logFile)
    
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    
    if(length(allCells) != 0){
      
      o <- getMatrixFromArrow(
        ArrowFile = ArrowFiles[x],
        useMatrix = useMatrix,
        useSeqnames = useSeqnames,
        cellNames = allCells, 
        ArchRProj = ArchRProj,
        verbose = FALSE,
        binarize = binarize,
        logFile = logFile
      )
      
      .logDiffTime(paste0("Completed ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
                   t1 = tstart, verbose = FALSE, logFile = logFile)
      
      o
      
    }else{
      
      NULL
      
    }
    
  }, threads = threads) 
  
  #ColData
  .logDiffTime("Organizing colData", t1 = tstart, verbose = verbose, logFile = logFile)
  cD <- lapply(seq_along(seL), function(x){
    colData(seL[[x]])
  }) %>% Reduce("rbind", .)
  
  #RowData
  .logDiffTime("Organizing rowData", t1 = tstart, verbose = verbose, logFile = logFile)
  rD1 <- rowData(seL[[1]])
  rD <- lapply(seq_along(seL), function(x){
    identical(rowData(seL[[x]]), rD1)
  }) %>% unlist %>% all
  if(!rD){
    stop("Error with rowData being equal for every sample!")
  }
  
  #RowRanges
  .logDiffTime("Organizing rowRanges", t1 = tstart, verbose = verbose, logFile = logFile)
  rR1 <- rowRanges(seL[[1]])
  rR <- lapply(seq_along(seL), function(x){
    identical(rowRanges(seL[[x]]), rR1)
  }) %>% unlist %>% all
  if(!rR){
    stop("Error with rowRanges being equal for every sample!")
  }
  
  #Assays
  nAssays <- names(assays(seL[[1]]))
  asy <- lapply(seq_along(nAssays), function(i){
    .logDiffTime(sprintf("Organizing Assays (%s of %s)", i, length(nAssays)), t1 = tstart, verbose = verbose, logFile = logFile)
    m <- lapply(seq_along(seL), function(j){
      assays(seL[[j]])[[nAssays[i]]]
    }) %>% Reduce("cbind", .)
    m
  }) %>% SimpleList()
  names(asy) <- nAssays
  
  .logDiffTime("Constructing SummarizedExperiment", t1 = tstart, verbose = verbose, logFile = logFile)
  if(!is.null(rR1)){
    se <- SummarizedExperiment(assays = asy, colData = cD, rowRanges = rR1)
    se <- sort(se)
  }else{
    se <- SummarizedExperiment(assays = asy, colData = cD, rowData = rD1)
  }
  rm(seL)
  gc()
  
  .logDiffTime("Finished Matrix Creation", t1 = tstart, verbose = verbose, logFile = logFile)
  
  se
  
}

#' Get a data matrix stored in an ArrowFile
#' 
#' This function gets a given data matrix from an individual ArrowFile.
#'
#' @param ArrowFile The path to an ArrowFile from which the selected data matrix should be obtained.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells from
#' the provided ArrowFile using `getCellNames()`.
#' @param ArchRProj An `ArchRProject` object to be used for getting additional information for cells in `cellColData`.
#' In some cases, data exists within the `ArchRProject` object that does not exist within the ArrowFiles. To access this data, you can
#' provide the `ArchRProject` object here.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromArrow <- function(
    ArrowFile = NULL, 
    useMatrix = "GeneScoreMatrix",
    useSeqnames = NULL,
    cellNames = NULL, 
    ArchRProj = NULL,
    verbose = TRUE,
    binarize = FALSE,
    logFile = createLogFile("getMatrixFromArrow")
){
  
  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = useMatrix, name = "useMatrix", valid = "character")
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromArrow Input-Parameters", logFile = logFile)
  
  ArrowFile <- .validArrow(ArrowFile)
  sampleName <- .sampleName(ArrowFile)
  
  seqnames <- .availableSeqnames(ArrowFile, subGroup = useMatrix)
  featureDF <- .getFeatureDF(ArrowFile, subGroup = useMatrix)
  .logThis(featureDF, paste0("featureDF ", sampleName), logFile = logFile)
  
  if(!is.null(useSeqnames)){
    seqnames <- seqnames[seqnames %in% useSeqnames]
  }
  
  if(length(seqnames) == 0){
    stop("No seqnames available!")
  }
  
  featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames), ]
  
  .logDiffTime(paste0("Getting ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
               t1 = tstart, verbose = verbose, logFile = logFile)
  
  if(!is.null(cellNames)){
    allCells <- .availableCells(ArrowFile = ArrowFile, subGroup = useMatrix)
    if(!all(cellNames %in% allCells)){
      stop("cellNames must all be within the ArrowFile!!!!")
    }
  }
  
  mat <- .getMatFromArrow(
    ArrowFile = ArrowFile, 
    featureDF = featureDF, 
    cellNames = cellNames, 
    useMatrix = useMatrix,
    binarize = binarize,
    useIndex = FALSE
  )
  .logThis(mat, paste0("mat ", sampleName), logFile = logFile)
  
  .logDiffTime(paste0("Organizing SE ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
               t1 = tstart, verbose = verbose, logFile = logFile)
  matrixClass <- h5read(ArrowFile, paste0(useMatrix, "/Info/Class"))
  
  if(matrixClass == "Sparse.Assays.Matrix"){
    rownames(mat) <- paste0(featureDF$name)
    splitIdx <- split(seq_len(nrow(mat)), featureDF$seqnames)
    mat <- lapply(seq_along(splitIdx), function(x){
      mat[splitIdx[[x]], , drop = FALSE]
    }) %>% SimpleList
    names(mat) <- names(splitIdx)
    featureDF <- featureDF[!duplicated(paste0(featureDF$name)), ,drop = FALSE]
    featureDF <- featureDF[,which(colnames(featureDF) %ni% "seqnames"), drop=FALSE]
    rownames(featureDF) <- paste0(featureDF$name)
  }else{
    mat <- SimpleList(mat)
    names(mat) <- useMatrix    
  }
  
  colData <- .getMetadata(ArrowFile)
  colData <- colData[colnames(mat[[1]]),,drop=FALSE]
  
  if(!is.null(ArchRProj)){
    projColData <- getCellColData(ArchRProj)[rownames(colData), ]
    colData <- cbind(colData, projColData[ ,colnames(projColData) %ni% colnames(colData)])
  }
  
  rowData <- tryCatch({
    makeGRangesFromDataFrame(featureDF, keep.extra.columns = TRUE)
  }, error = function(x){
    featureDF
  })
  
  se <- SummarizedExperiment(
    assays = mat,
    rowData = rowData,
    colData = colData
  )
  .logThis(se, paste0("se ", sampleName), logFile = logFile)
  
  se
  
}

.getMatFromArrow <- function(
    ArrowFile = NULL, 
    featureDF = NULL, 
    binarize = NULL, 
    cellNames = NULL,
    useMatrix = "TileMatrix", 
    useIndex = FALSE,
    threads = 1
){
  
  if(is.null(featureDF)){
    featureDF <- .getFeatureDF(ArrowFile, useMatrix)
  }
  
  if(any(c("seqnames","idx") %ni% colnames(featureDF))){
    stop("Need to provide featureDF with columns seqnames and idx!")
  }
  
  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  
  o <- h5closeAll()
  
  matClass <- h5read(ArrowFile, paste0(useMatrix,"/Info/Class"))
  if(matClass %ni% c("Sparse.Binary.Matrix", "Sparse.Integer.Matrix", "Sparse.Double.Matrix", "Sparse.Assays.Matrix")){
    stop("Arrow Mat is not a valid Sparse Matrix!")
  }
  if(is.null(binarize)){
    if(matClass == "Sparse.Binary.Matrix"){
      binarize <- TRUE
    }else{
      binarize <- FALSE
    }
  }
  if(matClass == "Sparse.Binary.Matrix"){
    if(!binarize){
      stop("Sparse Matrix in Arrow is Binarized! Set binarize = TRUE to use matrix!")
    }
  }
  
  matColNames <- paste0(.sampleName(ArrowFile), "#", h5read(ArrowFile, paste0(useMatrix,"/Info/CellNames")))
  if(!is.null(cellNames)){
    idxCols <- which(matColNames %in% cellNames)
  }else{
    idxCols <- seq_along(matColNames)
  }
  
  seqnames <- unique(featureDF$seqnames)
  
  mat <- .safelapply(seq_along(seqnames), function(x){
    
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex),]
    idxRows <- featureDFx$idx
    
    j <- Rle(
      values = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jValues")), 
      lengths = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jLengths"))
    )
    
    #Match J
    matchJ <- S4Vectors::match(j, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    if(useIndex){
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"), index = list(idxJ, 1))
    }else{
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"))[idxJ]
    }
    j <- matchJ[idxJ]
    
    #Match I
    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    i <- i[idxI]
    j <- j[idxI]
    i <- matchI[idxI]
    
    if(!binarize){
      x <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/x"))[idxJ][idxI]
    }else{
      x <- rep(1, length(j))
    }
    
    mat <- Matrix::sparseMatrix(
      i=as.vector(i),
      j=j,
      x=x,
      dims = c(length(idxRows), length(idxCols))
    )
    rownames(mat) <- rownames(featureDFx)
    
    rm(matchI, idxI, matchJ, idxJ, featureDFx, idxRows)
    
    return(mat)
    
  }, threads = threads) %>% Reduce("rbind", .)
  
  o <- h5closeAll()
  
  colnames(mat) <- matColNames[idxCols]
  
  #Double Check Order!
  mat <- mat[rownames(featureDF), , drop = FALSE]
  rownames(mat) <- NULL
  
  if(!is.null(cellNames)){
    mat <- mat[,cellNames,drop=FALSE]
  }
  
  return(mat)
  
}

####################################################################
# Helper read functioning
####################################################################
.getGroupMatrix <- function(
    ArrowFiles = NULL, 
    featureDF = NULL, 
    groupList = NULL,
    threads = 1, 
    useIndex = FALSE, 
    verbose = TRUE, 
    useMatrix = "TileMatrix",
    asSparse = FALSE,
    tstart = NULL
){
  
  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  #########################################
  # Construct Matrix
  #########################################
  seqnames <- unique(featureDF$seqnames)
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  cellNames <- unlist(groupList, use.names = FALSE) ### UNIQUE here? doublet check QQQ
  
  allCellsList <- lapply(seq_along(ArrowFiles), function(x){
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    if(length(allCells) != 0){
      allCells
    }else{
      NULL
    }
  })
  
  mat <- .safelapply(seq_along(seqnames), function(x){
    
    .logDiffTime(sprintf("Constructing Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)
    
    #Construct Matrix
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex), ]
    
    matChr <- matrix(0, nrow = nrow(featureDFx), ncol = length(groupList))
    colnames(matChr) <- names(groupList)
    rownames(matChr) <- rownames(featureDFx)
    
    for(y in seq_along(ArrowFiles)){
      
      allCells <- allCellsList[[y]]
      
      if(!is.null(allCells)){
        
        maty <- .getMatFromArrow(
          ArrowFile = ArrowFiles[y], 
          useMatrix = useMatrix,
          featureDF = featureDFx, 
          cellNames = allCells, 
          useIndex = useIndex
        )
        
        for(z in seq_along(groupList)){
          
          #Check Cells In Group
          cellsGroupz <- groupList[[z]]
          idx <- BiocGenerics::which(colnames(maty) %in% cellsGroupz)
          
          #If In Group RowSums
          if(length(idx) > 0){
            matChr[,z] <- matChr[,z] + Matrix::rowSums(maty[,idx,drop=FALSE])
          }
          
        }
        
        rm(maty)
        
      }
      
      
      if(y %% 20 == 0 | y %% length(ArrowFiles) == 0){
        gc()
      } 
      
    }
    
    if(asSparse){
      matChr <- as(matChr, "dgCMatrix")
    }
    
    .logDiffTime(sprintf("Finished Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)
    
    matChr
    
  }, threads = threads) %>% Reduce("rbind", .)
  
  mat <- mat[rownames(featureDF), , drop = FALSE]
  
  .logDiffTime("Successfully Created Group Matrix", tstart, verbose = verbose)
  
  gc()
  
  return(mat)
  
}

.getPartialMatrix <- function(
    ArrowFiles = NULL, 
    featureDF = NULL, 
    cellNames = NULL, 
    progress = TRUE, 
    threads = 1, 
    useMatrix = "TileMatrix",
    doSampleCells = FALSE, 
    sampledCellNames = NULL, 
    tmpPath = .tempfile(pattern = paste0("tmp-partial-mat")), 
    useIndex = FALSE,
    tstart = NULL,
    verbose = TRUE,
    logFile = NULL
){
  
  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  #########################################
  # Construct Matrix
  #########################################
  
  mat <- .safelapply(seq_along(ArrowFiles), function(x){
    
    .logDiffTime(sprintf("Getting Partial Matrix %s of %s", x, length(ArrowFiles)), tstart, verbose = verbose)
    
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    
    if(length(allCells) == 0){
      if(doSampleCells){
        return(list(mat = NULL, out = NULL))
      }else{
        return(NULL)
      }
    }
    
    o <- h5closeAll()
    matx <- .getMatFromArrow(
      ArrowFile = ArrowFiles[x], 
      featureDF = featureDF, 
      cellNames = allCells,
      useMatrix = useMatrix, 
      useIndex = useIndex
    )
    
    if(doSampleCells){
      
      #Save Temporary Matrix
      outx <- paste0(tmpPath, "-", .sampleName(ArrowFiles[x]), ".rds")
      .safeSaveRDS(matx, outx, compress = FALSE)     
      
      #Sample Matrix 
      matx <- matx[, which(colnames(matx) %in% sampledCellNames),drop = FALSE]
      
      return(list(mat = matx, out = outx))
      
    }else{
      
      return(matx)
      
    }
    
  }, threads = threads)
  
  gc()
  
  
  if(doSampleCells){
    
    matFiles <- lapply(mat, function(x) x[[2]]) %>% Reduce("c", .)
    mat <- lapply(mat, function(x) x[[1]]) %>% Reduce("cbind", .)
    if(!all(sampledCellNames %in% colnames(mat))){
      .logThis(sampledCellNames, "cellNames supplied", logFile = logFile)
      .logThis(colnames(mat), "cellNames from matrix", logFile = logFile)
      stop("Error not all cellNames found in partialMatrix")
    }
    mat <- mat[,sampledCellNames, drop = FALSE]
    mat <- .checkSparseMatrix(mat, length(sampledCellNames))
    
    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)
    
    return(list(mat = mat, matFiles = matFiles))
    
  }else{
    
    mat <- Reduce("cbind", mat)
    if(!all(cellNames %in% colnames(mat))){
      .logThis(cellNames, "cellNames supplied", logFile = logFile)
      .logThis(colnames(mat), "cellNames from matrix", logFile = logFile)
      stop("Error not all cellNames found in partialMatrix")
    }
    mat <- mat[,cellNames, drop = FALSE]
    mat <- .checkSparseMatrix(mat, length(cellNames))
    
    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)
    
    return(mat)
    
  }
  
  
}

.checkSparseMatrix <- function(x, ncol = NULL){
  isSM <- is(x, 'sparseMatrix')
  if(!isSM){
    if(is.null(ncol)){
      stop("ncol must not be NULL if x is not a matrix!")
    }
    cnames <- tryCatch({
      names(x)
    }, error = function(e){
      colnames(x)
    })
    if(length(cnames) != ncol){
      stop("cnames != ncol!")
    }
    x <- Matrix::Matrix(matrix(x, ncol = ncol), sparse=TRUE)
    colnames(x) <- cnames
  }
  x
}

########################################################################
# Compute Summary Statistics!
########################################################################

.getRowSums <- function(
    ArrowFiles = NULL,
    useMatrix = NULL,
    seqnames = NULL,
    verbose = TRUE,
    tstart = NULL,
    filter0 = FALSE,
    threads = 1,
    addInfo = FALSE
){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  if(is.null(seqnames)){
    seqnames <- .availableSeqnames(ArrowFiles, useMatrix)
  }
  
  #Compute RowSums
  summaryDF <- .safelapply(seq_along(seqnames), function(x){
    o <- h5closeAll()
    for(y in seq_along(ArrowFiles)){
      if(y == 1){
        sumy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
      }else{
        sumy1 <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
        if(length(sumy1) != length(sumy)){
          stop("rowSums lengths do not match in ArrowFiles for a seqname!")
        }else{
          sumy <- sumy + sumy1
        }
      }
    }
    #Return Setup In Feature DF Format (seqnames, idx columns)
    DataFrame(seqnames = Rle(seqnames[x], lengths = length(sumy)), idx = seq_along(sumy), rowSums = as.vector(sumy))
  }, threads = threads) %>% Reduce("rbind", .)
  
  if(addInfo){
    featureDF <- .getFeatureDF(ArrowFiles, useMatrix)
    rownames(featureDF) <- paste0(featureDF$seqnames, "_", featureDF$idx)
    rownames(summaryDF) <- paste0(summaryDF$seqnames, "_", summaryDF$idx)
    featureDF <- featureDF[rownames(summaryDF), , drop = FALSE]
    featureDF$rowSums <- summaryDF[rownames(featureDF), "rowSums"]
    summaryDF <- featureDF
    rownames(summaryDF) <- NULL
    remove(featureDF)
  }
  
  if(filter0){
    summaryDF <- summaryDF[which(summaryDF$rowSums > 0), ,drop = FALSE]
  }
  
  return(summaryDF)
  
}

.getRowVars <- function(
    ArrowFiles = NULL,
    seqnames = NULL,
    useMatrix = NULL,
    useLog2 = FALSE,
    threads = 1
){
  
  .combineVariances <- function(dfMeans = NULL, dfVars = NULL, ns = NULL){
    
    #https://rdrr.io/cran/fishmethods/src/R/combinevar.R
    
    if(ncol(dfMeans) != ncol(dfVars) | ncol(dfMeans) != length(ns)){
      stop("Means Variances and Ns lengths not identical")
    }
    
    #Check if samples have NAs due to N = 1 sample or some other weird thing.
    #Set it to min non NA variance
    dfVars <- lapply(seq_len(nrow(dfVars)), function(x){
      vx <- dfVars[x, , drop = FALSE]
      if(any(is.na(vx))){
        vx[is.na(vx)] <- min(vx[!is.na(vx)])
      }
      vx
    }) %>% Reduce("rbind", .)
    
    combinedMeans <- rowSums(t(t(dfMeans) * ns)) / sum(ns)
    summedVars <- rowSums(t(t(dfVars) * (ns - 1)) + t(t(dfMeans^2) * ns))
    combinedVars <- (summedVars - sum(ns)*combinedMeans^2)/(sum(ns)-1)
    
    data.frame(combinedVars = combinedVars, combinedMeans = combinedMeans)
    
  }
  
  featureDF <- .getFeatureDF(ArrowFiles, useMatrix)
  
  if(!is.null(seqnames)){
    featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames),]
  }
  
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  fnames <- rownames(featureDF)
  
  featureDF <- S4Vectors::split(featureDF, as.character(featureDF$seqnames))
  
  ns <- lapply(seq_along(ArrowFiles), function(y){
    length(.availableCells(ArrowFiles[y], useMatrix))
  }) %>% unlist
  
  #Compute RowVars
  summaryDF <- .safelapply(seq_along(featureDF), function(x){
    
    o <- h5closeAll()
    seqx <- names(featureDF)[x]
    meanx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    varx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    
    for(y in seq_along(ArrowFiles)){
      
      if(useLog2){
        meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeansLog2"))
        varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVarsLog2")) 
      }else{
        meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeans"))
        varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVars"))
      }
      
    }
    
    cbind(featureDF[[x]], DataFrame(.combineVariances(meanx, varx, ns)))
    
  }, threads = threads) %>% Reduce("rbind", .)
  
  summaryDF <- summaryDF[fnames, , drop = FALSE]
  
  return(summaryDF)
  
}

.getColSums <- function(
    ArrowFiles = NULL,
    seqnames = NULL,
    useMatrix = NULL,
    verbose = TRUE,
    tstart = NULL,
    threads = 1
){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  #Compute ColSums
  cS <- .safelapply(seq_along(seqnames), function(x){
    
    lapply(seq_along(ArrowFiles), function(y){
      
      o <- h5closeAll()
      cSy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/colSums"))
      rownames(cSy) <- .availableCells(ArrowFiles[y], useMatrix)
      cSy
      
    }) %>% Reduce("rbind", .)
    
  }, threads = threads) %>% Reduce("cbind", .) %>% rowSums
  
  .logDiffTime("Successfully Computed colSums", tstart, verbose = verbose)
  
  return(cS)
  
}

# h5read implementation for optimal reading
.h5read <- function(
    file = NULL,
    name = NULL,
    method = "fast",
    index = NULL,
    start = NULL,
    block = NULL,
    count = NULL
){
  
  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- H5Fopen(file, "H5F_ACC_RDONLY")
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}

.getMatrixClass <- function(
    ArrowFiles = NULL, 
    useMatrix = NULL,
    threads = getArchRThreads()
){
  
  threads <- min(length(ArrowFiles), threads)
  
  matrixClass <- .safelapply(seq_along(ArrowFiles), function(i){
    h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Class"))
  }, threads = threads) %>% unlist %>% unique
  
  if(length(matrixClass) != 1){
    stop("Not all matrix classes are the same!")
  }
  
  matrixClass
  
}

.getMatrixUnits <- function(
    ArrowFiles = NULL, 
    useMatrix = NULL,
    threads = getArchRThreads()
){
  
  threads <- min(length(ArrowFiles), threads)
  
  matrixUnits <- .safelapply(seq_along(ArrowFiles), function(i){
    tryCatch({ #This handles backwards compatibility!
      h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Units"))
    }, error = function(x){
      "None"
    })
  }, threads = threads) %>% unlist %>% unique
  
  if(length(matrixUnits) != 1){
    stop("Not all matrix units are the same!")
  }
  
  matrixUnits
  
}


.initializeMat <- function(
    ArrowFile = NULL,
    Group = NULL,
    Class = "Double",
    Units = "none",
    cellNames = NULL,
    featureDF = NULL,
    params = NULL,
    date = Sys.Date(),
    force = FALSE
){
  
  #Add Group Entry of SparseMatrix Format
  #This Includes the following format
  #
  # Info 
  #   - Class - Sparse.Integer.Matrix = Sparse Matrix with Integer Entries
  #           - Sparse.Binary.Matrix = Sparse Matrix with Binary ie no x values
  #           - Sparse.Double.Matrix = Sparse Matrix with Double/Numeric Entries
  #           - Sparse.Assays.Matrix = Sparse Matrices with same feature names (same cell x feature)
  #                                    separated by different seqnames (ie assayNames)
  #   - CellNames ie Colnames
  #   - FeatureDF dataframe that describes the rows of each seqname
  #   - Params Params that are used for construction to be checked when comparing Arrows
  #   - Date Date of Creation
  # Chr1
  #   - i, j (as an Rle), x, and rowSums,colSums,rowVars,etc.
  # Chr2
  # Chr3
  # ..
  #
  
  if(!suppressMessages(h5createGroup(ArrowFile, paste0(Group)))){
    if(force){
      h5delete(ArrowFile, paste0(Group))
      h5createGroup(ArrowFile, paste0(Group))
    }else{
      stop("Matrix Group Already Exists! Set force = TRUE to overwrite!")
    }
  }
  o <- h5createGroup(ArrowFile, paste0(Group, "/Info"))
  
  if(tolower(Class)=="binary"){
    
    o <- h5write(obj = "Sparse.Binary.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))
    
  }else if(tolower(Class)=="integer"){
    
    o <- h5write(obj = "Sparse.Integer.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))
    
  }else if(tolower(Class)=="double"){
    
    o <- h5write(obj = "Sparse.Double.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))
    
  }else if(tolower(Class)=="assays"){
    
    o <- h5write(obj = "Sparse.Assays.Matrix", file = ArrowFile, name = paste0(Group, "/Info/Class"))
    
  }else{
    
    stop("Matrix Class Not Supported!")
    
  }
  
  ##########
  # Add Units To Class
  ##########
  if(!is.character(Units)){
    stop("Please provide Units as character when writing matrix to Arrow!")
  }
  o <- h5write(obj = Units, file = ArrowFile, name = paste0(Group, "/Info/Units"))  
  
  ##########
  # Cell Names in Arrow
  ##########
  splitNames <- stringr::str_split(cellNames, pattern = "#", simplify=TRUE)
  if(ncol(splitNames) > 2){
    stop("Found error with cell names containing multiple # characters!")
  }else{
    cellNames <- splitNames[,ncol(splitNames)]
  }
  o <- h5write(obj = cellNames, file = ArrowFile, name = paste0(Group,"/Info/CellNames"))
  
  ##########
  # FeatureDF in Arrow
  ##########
  df <- data.frame(featureDF, stringsAsFactors = FALSE)
  stopifnot(all(c("seqnames","idx") %in% colnames(featureDF)))
  o <- h5write(obj = df, file = ArrowFile, name = paste0(Group,"/Info/FeatureDF"))
  
  ##########
  # Parameters for Matrix for Validity in Arrow
  ##########  
  o <- h5write(obj = params, file = ArrowFile, name = paste0(Group,"/Info/Params"))
  
  ##########
  # Date of Creation
  ##########  
  o <- h5write(obj = paste0(date), file = ArrowFile, name = paste0(Group,"/Info/Date"))
  
  return(0)
  
}

.addMatToArrow <- function(
    mat = NULL, 
    ArrowFile = NULL,
    Group = NULL,
    binarize = FALSE, 
    addRowSums = FALSE,
    addColSums = FALSE,
    addRowMeans = FALSE,
    addRowVars = FALSE,
    addRowVarsLog2 = FALSE,
    logFile = NULL
){
  
  stopifnot(inherits(mat, "dgCMatrix"))
  
  checkCells <- .availableCells(ArrowFile, dirname(Group))
  if(!identical(paste0(colnames(mat)), paste0(checkCells))){
    .logThis(colnames(mat), "colnames", logFile=logFile)
    .logThis(checkCells, "cellNames", logFile=logFile)
    .logMessage(paste0("Mismatch = ", sum(colnames(mat) != checkCells)))
    .logMessage("CellNames in Matrix Group do not Match CellNames in Matrix Being Written!",logFile=logFile)
    stop("CellNames in Matrix Group do not Match CellNames in Matrix Being Written!")
  }
  
  #Create Group
  o <- h5closeAll()
  o <- h5createGroup(ArrowFile, Group)
  
  #Convert Columns to Rle
  j <- Rle(findInterval(seq(mat@x)-1,mat@p[-1]) + 1)
  
  #Info
  lengthRle <- length(j@lengths)
  lengthI <- length(mat@i)
  
  #Create Data Set
  o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group,"/i"), storage.mode = "integer", 
                                    dims = c(lengthI, 1), level = 0))
  
  o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group,"/jLengths"), storage.mode = "integer", 
                                    dims = c(lengthRle, 1), level = 0))
  
  o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group,"/jValues"), storage.mode = "integer", 
                                    dims = c(lengthRle, 1), level = 0))
  
  #Write Data Set
  o <- .suppressAll(h5write(obj = mat@i + 1, file = ArrowFile, name = paste0(Group,"/i")))
  o <- .suppressAll(h5write(obj = j@lengths, file = ArrowFile, name = paste0(Group,"/jLengths")))
  o <- .suppressAll(h5write(obj = j@values, file = ArrowFile, name = paste0(Group,"/jValues")))
  
  #If binary dont store x
  if(!binarize){
    
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/x"), storage.mode = "double", 
                                      dims = c(lengthI, 1), level = 0))
    
    o <- .suppressAll(h5write(obj = mat@x, file = ArrowFile, name = paste0(Group, "/x")))
    
  }else{
    
    mat@x[mat@x > 0] <- 1
    
  }
  
  if(addColSums){
    cS <- Matrix::colSums(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/colSums"), storage.mode = "double", 
                                      dims = c(ncol(mat), 1), level = 0))
    o <- .suppressAll(h5write(obj = cS, file = ArrowFile, name = paste0(Group, "/colSums")))
    
  }
  
  if(addRowSums){
    rS <- Matrix::rowSums(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowSums"), storage.mode = "double", 
                                      dims = c(nrow(mat), 1), level = 0))
    o <- .suppressAll(h5write(obj = rS, file = ArrowFile, name = paste0(Group, "/rowSums")))
    
  }
  
  if(addRowMeans){
    rM <- Matrix::rowMeans(mat)
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeans"), storage.mode = "double", 
                                      dims = c(nrow(mat), 1), level = 0))
    o <- .suppressAll(h5write(obj = rM, file = ArrowFile, name = paste0(Group, "/rowMeans")))
    
  }
  
  if(addRowVars){
    if(!addRowMeans){
      rM <- Matrix::rowMeans(mat)
    }
    rV <- computeSparseRowVariances(mat@i + 1, mat@x, rM, n = ncol(mat))
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVars"), storage.mode = "double", 
                                      dims = c(nrow(mat), 1), level = 0))
    o <- .suppressAll(h5write(obj = rV, file = ArrowFile, name = paste0(Group, "/rowVars")))
    
  }
  
  if(addRowVarsLog2){
    
    mat@x <- log2(mat@x + 1) #log-normalize
    rM <- Matrix::rowMeans(mat)
    idx <- seq_along(rM)
    idxSplit <- .splitEvery(idx, 2000)
    
    #Make sure not too much memory so split into 2000 gene chunks
    #Check this doesnt cause memory mapping issues!
    rV <- lapply(seq_along(idxSplit), function(x){
      idxX <- idxSplit[[x]]
      matx <- mat[idxX, , drop = FALSE]
      computeSparseRowVariances(matx@i + 1, matx@x, rM[idxX], n = ncol(mat))
    }) %>% unlist
    
    #Have to write rowMeansLog2 as well
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeansLog2"), storage.mode = "double", 
                                      dims = c(nrow(mat), 1), level = 0))
    o <- .suppressAll(h5write(obj = rM, file = ArrowFile, name = paste0(Group, "/rowMeansLog2")))
    
    #Write rowVarsLog2
    o <- .suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVarsLog2"), storage.mode = "double", 
                                      dims = c(nrow(mat), 1), level = 0))
    o <- .suppressAll(h5write(obj = rV, file = ArrowFile, name = paste0(Group, "/rowVarsLog2")))
  }
  
  #Clean Up Memorys
  rm(j,mat)
  gc()
  
  o <- h5closeAll()
  
  return(0)
  
}