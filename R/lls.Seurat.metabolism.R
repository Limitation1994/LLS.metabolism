#' @title Metabolic Pathway Activity Scoring
#' @description This function performs metabolic pathway activity scoring using various methods.
#' @param sce A Seurat object containing RNA assay data.
#' @param method The method to use for scoring ("VISION", "AUCell", "ssGSEA", "gsva").
#' @param missing_values A boolean indicating whether to fill missing values.
#' @param ncores Number of cores to use for parallel processing.
#' @param metabolism.type Type of metabolism pathway to use ("KEGG", "REACTOME", "LLS").
#' @return A Seurat object with metabolic pathway scores added.
#' @export
lls.Seurat.metabolism <- function(sce, method = "VISION", missing_values = F, ncores = 2, metabolism.type = "LLS") {
  library(VISION)
  library(AUCell)
  library(GSEABase)
  library(GSVA)
  library(progress)

  countexp <- sce@assays$RNA@counts
  countexp <- data.frame(as.matrix(countexp))

  metabolism_LLS_genesets <- system.file("data", "metabolism_LLS.gmt", package = "LLS.metabolism")
  metabolism_KEGG_genesets <- system.file("data", "metabolism_KEGG.gmt", package = "LLS.metabolism")
  metabolism_REACTOME_genesets <- system.file("data", "metabolism_REACTOME.gmt", package = "LLS.metabolism")

  if (metabolism.type == "KEGG") {
    gmtFile <- metabolism_KEGG_genesets
    cat("KEGG metabolic pathway will be used for subsequent analysis\n")
  }
  if (metabolism.type == "REACTOME") {
    gmtFile <- metabolism_REACTOME_genesets
    cat("REACTOME metabolic pathway will be used for subsequent analysis\n")
  }
  if (metabolism.type == "LLS") {
    gmtFile <- metabolism_LLS_genesets
    cat("LLS metabolic pathway will be used for subsequent analysis\n")
  }

  if (!missing_values) {
    count_exp <- countexp
  } else {
    cat("Filling missing values\n")
    result.completed <- alra(as.matrix(countexp))
    count_exp <- result.completed[[3]]
    row.names(count_exp) <- row.names(countexp)
  }

  cat("Metabolic pathway activity scoring in progress\n")

  # VISION
  if (method == "VISION") {
    cat("VISION-----Metabolic pathway activity scoring in progress\n")

    pb <- progress_bar$new(
      format = "  [:bar] :percent :elapsedfull",
      total = 5,  # 总的步骤数
      clear = FALSE
    )

    n.umi <- colSums(count_exp)
    pb$tick()
    cat("Step 1: Calculated total UMI counts.\n")

    scaled_counts <- t(t(count_exp) / n.umi) * median(n.umi)
    pb$tick()
    cat("Step 2: Scaled count matrix based on UMI counts.\n")

    vis <- Vision(scaled_counts, signatures = gmtFile)
    pb$tick()
    cat("Step 3: Created Vision object with scaled counts and signatures.\n")

    options(mc.cores = ncores)
    vis <- analyze(vis)
    pb$tick()
    cat("Step 4: Analyzed Vision object with multiple cores---Slightly slow, please wait patiently.\n")

    signature_exp <- data.frame(t(vis@SigScores))
    pb$tick()
    cat("Step 5: Transformed signature scores into data frame.\n")
  }

  # AUCell
  if (method == "AUCell") {
    cat("AUCell-----Pathway activity scoring in progress\n")

    pb <- progress_bar$new(
      format = "  [:bar] :percent :elapsedfull",
      total = 4,  # 总的步骤数
      clear = FALSE
    )

    cells_rankings <- AUCell_buildRankings(as.matrix(count_exp), nCores = ncores, plotStats = F)
    pb$tick()
    cat("Step 1: Built cell rankings.\n")

    geneSets <- getGmt(gmtFile)
    pb$tick()
    cat("Step 2: Loaded gene sets from GMT file.\n")

    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    pb$tick()
    cat("Step 3: Calculated AUC scores for cells.\n")

    signature_exp <- data.frame(getAUC(cells_AUC))
    pb$tick()
    cat("Step 4: Transformed AUC scores into data frame.\n")
  }

  # ssGSEA
  if (method == "ssGSEA") {
    cat("ssGSEA-----Pathway activity scoring in progress\n")

    pb <- progress_bar$new(
      format = "  [:bar] :percent :elapsedfull",
      total = 4,  # 总的步骤数
      clear = FALSE
    )

    geneSets <- getGmt(gmtFile)
    pb$tick()
    cat("Step 1: Loaded gene sets from GMT file.\n")

    gsva_es <- gsva(as.matrix(count_exp), geneSets, method = "ssgsea", kcdf = "Poisson", parallel.sz = ncores)
    pb$tick()
    cat("Step 2: Calculated ssGSEA enrichment scores.\n")

    signature_exp <- data.frame(gsva_es)
    pb$tick()
    cat("Step 3: Transformed enrichment scores into data frame.\n")
  }

  # GSVA
  if (method == "gsva") {
    cat("GSVA-----Pathway activity scoring in progress\n")

    pb <- progress_bar$new(
      format = "  [:bar] :percent :elapsedfull",
      total = 4,  # 总的步骤数
      clear = FALSE
    )

    geneSets <- getGmt(gmtFile)
    pb$tick()
    cat("Step 1: Loaded gene sets from GMT file.\n")

    gsva_es <- gsva(as.matrix(count_exp), geneSets, method = "gsva", kcdf = "Poisson", parallel.sz = ncores)
    pb$tick()
    cat("Step 2: Calculated GSVA enrichment scores.\n")

    signature_exp <- data.frame(gsva_es)
    pb$tick()
    cat("Step 3: Transformed enrichment scores into data frame.\n")
  }

  cat("Thank you very much for using our R package and hope it is helpful to you. If possible, the relevant citations are as follows: AAA\n")

  sce@assays$METABOLISM$score <- signature_exp
  sce
}
