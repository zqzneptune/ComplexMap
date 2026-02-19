utils::globalVariables(c("domainId", "Description", "p.adjust", "complexId",
                         "primaryFunctionalDomain", "topEnrichedFunctions",
                         "colorHex", "proteinCount", "Fold", "score", "as",
                         "domainLabel", "topIdx", "baseColor", "head"))

#' Generate Node Attributes for a Complex Network
#'
#' @description
#' Creates a detailed attribute table for each complex, emphasizing functional
#' specificity for network visualization.
#'
#' @details
#' This function performs several steps to generate rich node attributes:
#' 1.  It aggregates all enriched terms from the input `enrichments` list.
#' 2.  It calculates a **Specificity Score** for each term using the formula:
#'     `-log10(p.adjust) * log2(Fold)`. This ensures that specific, high-fold
#'     enrichment terms are prioritized over broad, generic terms.
#' 3.  A term-complex matrix is built, and terms are clustered. The clustering
#'     uses **Average Linkage** and forces a higher number of clusters to
#'     preserve functional diversity (avoiding "monochromatic" maps).
#' 4.  A unique color is assigned to each functional domain.
#' 5.  For each complex, the **Primary Functional Domain** is assigned to the
#'     enriched term with the highest Specificity Score.
#' 6.  A unique "blended" color is calculated by mixing domain colors, weighted
#'     by the Specificity Score.
#'
#' @param complexes A named list of protein complexes.
#' @param enrichments A named list of enrichment results. Must contain 'Fold'
#'   column for specificity weighting.
#' @param geneSetDb Optional named list of gene sets for semantic clustering.
#' @param similarityMethod Distance method for clustering ("jaccard", "overlap", etc.).
#'   Defaults to "jaccard" to penalize size differences.
#' @param verbose Logical.
#'
#' @return A `tibble` with complex attributes.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
generateNodeAttributes <- function(complexes, enrichments,
                                   geneSetDb = NULL,
                                   similarityMethod = "jaccard",
                                   verbose = TRUE) {
  if (verbose) {
    message("Generating node attributes (prioritizing functional specificity)...")
  }
  
  allTermsDf <- dplyr::bind_rows(enrichments, .id = "complexId")
  
  if (nrow(allTermsDf) == 0) {
    if (verbose) message("No enriched terms found; returning basic attributes.")
    proteinStr <- vapply(complexes, paste, collapse = ",", FUN.VALUE = character(1))
    return(tibble::tibble(
      complexId = names(complexes),
      proteinCount = lengths(complexes),
      proteins = proteinStr,
      primaryFunctionalDomain = "Unenriched",
      topEnrichedFunctions = NA_character_,
      colorHex = "#CCCCCC"
    ))
  }
  
  enrichedComplexIds <- names(enrichments)
  allTerms <- unique(allTermsDf$ID)
  
  # --- 1. Cluster Terms (Force Diversity) ---
  if (!is.null(geneSetDb)) {
    if (verbose) message(sprintf("    -> Clustering %d terms using gene set similarity (%s)", length(allTerms), similarityMethod))
    distMatrix <- .computeGeneSetDistance(allTerms, geneSetDb, similarityMethod)
  } else {
    if (verbose) message(sprintf("    -> Clustering %d terms using co-occurrence (%s)", length(allTerms), similarityMethod))
    distMatrix <- .computeCooccurrenceDistance(allTerms, enrichments, enrichedComplexIds, similarityMethod)
  }
  
  # Use the REFACTORED .clusterTermsOptimal
  clusterResult <- .clusterTermsOptimal(distMatrix, allTerms, verbose)
  termDomains <- clusterResult$termDomains
  numDomains <- clusterResult$numDomains
  
  termToDomainMap <- tibble::tibble(term = names(termDomains), domainId = termDomains)
  
  # --- 2. Create Domain Labels ---
  # Pick the label for the Domain (for the Legend) based on specificity (Fold)
  domainLabels <- allTermsDf %>%
    dplyr::left_join(termToDomainMap, by = c("ID" = "term")) %>%
    dplyr::group_by(domainId) %>%
    dplyr::summarise(
      # Pick the most "specific" description (highest Fold) to label the color in the legend
      domainLabel = Description[which.max(Fold)][1], 
      .groups = "drop"
    )
  
  domainColors <- .assignDistinctColors(numDomains) %>%
    dplyr::mutate(domainId = seq_len(numDomains)) %>%
    dplyr::left_join(domainLabels, by = "domainId")
  
  # --- 3. Calculate Per-Complex Attributes (Weighted by Specificity) ---
  # [REFACTOR START] -----------------------------------------------------------
  # Logic changed: Instead of blending colors, we now strictly assign the color 
  # of the dominant functional cluster. This ensures 1:1 mapping between 
  # Primary Function and Color.
  
  allEnrichDf <- dplyr::bind_rows(lapply(names(enrichments), function(cid) {
    df <- enrichments[[cid]]
    if (!"p.adjust" %in% names(df)) df$p.adjust <- NA_real_
    if (!"Fold" %in% names(df)) df$Fold <- 1.0
    
    # Calculate Specificity Score
    df %>% dplyr::mutate(
      complexId = cid,
      score = -log10(pmax(p.adjust, 1e-300)) * log2(pmax(Fold, 1.1))
    )
  })) %>%
    dplyr::left_join(termToDomainMap, by = c("ID" = "term")) %>%
    # Join domainColors here so every term knows its Cluster Label and Base Color
    dplyr::left_join(domainColors, by = "domainId")
  
  complexSummary <- allEnrichDf %>%
    dplyr::group_by(complexId) %>%
    dplyr::summarise(
      # Find the index of the highest scoring term for this complex
      topIdx = which.max(score)[1],
      
      # Assign Primary Function: Use the DOMAIN LABEL (Cluster Representative).
      # This ensures that even if the specific terms differ slightly, if they 
      # belong to the same cluster, they get the same Label and Color.
      primaryFunctionalDomain = domainLabel[topIdx],
      
      # Assign Color: Strictly use the Base Color of that Domain.
      # No blending. This guarantees consistency with the primary function.
      colorHex = baseColor[topIdx],
      
      # Keep specific details in the 'topEnrichedFunctions' text
      topEnrichedFunctions = paste(unique(head(Description[order(score, decreasing = TRUE)], 3)), collapse = "; "),
      
      .groups = "drop"
    )
  # [REFACTOR END] -------------------------------------------------------------
  
  
  # --- 4. Finalize ---
  complexMeta <- tibble::tibble(
    complexId = names(complexes),
    proteinCount = lengths(complexes),
    proteins = vapply(complexes, paste, collapse = ",", FUN.VALUE = character(1))
  )
  
  finalDf <- complexMeta %>%
    dplyr::left_join(complexSummary, by = "complexId") %>%
    dplyr::mutate(
      primaryFunctionalDomain = dplyr::coalesce(primaryFunctionalDomain, "Unenriched"),
      topEnrichedFunctions = dplyr::coalesce(topEnrichedFunctions, NA_character_),
      colorHex = dplyr::coalesce(colorHex, "#CCCCCC")
    )
  
  return(finalDf)
}

# === INTERNAL HELPER FUNCTIONS ===

#' @keywords internal
.computeGeneSetDistance <- function(terms, geneSetDb, method) {
  validTerms <- intersect(terms, names(geneSetDb))
  if (length(validTerms) == 0) stop("No terms found in gene set database")
  
  termGeneSets <- geneSetDb[validTerms]
  nTerms <- length(termGeneSets)
  if (nTerms == 1) return(stats::as.dist(matrix(0, 1, 1, dimnames = list(validTerms, validTerms))))
  
  allGenes <- unique(unlist(termGeneSets, use.names = FALSE))
  # Sparse matrix approach for speed
  geneTermMatrix <- Matrix::sparseMatrix(
    i = match(unlist(termGeneSets), allGenes),
    j = rep(seq_along(termGeneSets), lengths(termGeneSets)),
    x = 1, dims = c(length(allGenes), nTerms),
    dimnames = list(allGenes, validTerms)
  )
  
  intersectionMatrix <- as.matrix(Matrix::crossprod(geneTermMatrix)) # T(M) %*% M
  setSizes <- Matrix::colSums(geneTermMatrix)
  
  simMatrix <- matrix(0, nrow = nTerms, ncol = nTerms, dimnames = list(validTerms, validTerms))
  
  # Vectorized Jaccard/Overlap/Dice
  if (method == "jaccard") {
    unionMatrix <- outer(setSizes, setSizes, "+") - intersectionMatrix
    simMatrix <- intersectionMatrix / unionMatrix
  } else if (method == "overlap") {
    minMatrix <- outer(setSizes, setSizes, pmin)
    simMatrix <- intersectionMatrix / minMatrix
  } else if (method == "cosine") {
    sqrtMatrix <- outer(setSizes, setSizes, "*")
    simMatrix <- intersectionMatrix / sqrt(sqrtMatrix)
  } else {
    # Default to Jaccard if unknown or "dice" if requested manually
    if (method == "dice") {
      simMatrix <- (2 * intersectionMatrix) / outer(setSizes, setSizes, "+")
    } else {
      unionMatrix <- outer(setSizes, setSizes, "+") - intersectionMatrix
      simMatrix <- intersectionMatrix / unionMatrix
    }
  }
  
  diag(simMatrix) <- 1
  stats::as.dist(1 - simMatrix)
}

#' @keywords internal
.computeCooccurrenceDistance <- function(terms, enrichments, complexIds, method) {
  termComplexMatrix <- matrix(0L, nrow = length(terms), ncol = length(complexIds),
                              dimnames = list(terms, complexIds))
  idsList <- lapply(enrichments, function(x) x$ID)
  termIdx <- match(unlist(idsList), terms)
  complexIdx <- rep(seq_along(complexIds), lengths(idsList))
  termComplexMatrix[cbind(termIdx, complexIdx)] <- 1L
  
  if (nrow(termComplexMatrix) <= 1) return(stats::as.dist(0))
  
  if (requireNamespace("philentropy", quietly = TRUE) && 
      method %in% c("jaccard", "overlap", "cosine", "dice")) {
    # FIX: Wrap philentropy matrix in as.dist()
    raw_dist <- philentropy::distance(termComplexMatrix, method = method, use.row.names = TRUE)
    return(stats::as.dist(raw_dist))
  }
  
  # Fallback to binary distance if philentropy is missing or method unknown
  return(stats::dist(termComplexMatrix, method = "binary"))
}

#' @keywords internal
.clusterTermsOptimal <- function(distMatrix, terms, verbose) {
  nTerms <- length(terms)
  if (nTerms <= 1) {
    return(list(termDomains = stats::setNames(1, terms), numDomains = 1))
  }
  
  # --- REFACTORED LOGIC FOR DIVERSITY ---
  # 1. Use Average Linkage to preserve outliers (branches) and avoid "blobbing"
  # This requires distMatrix to be a 'dist' object (FIXED in computeCooccurrenceDistance)
  hc <- stats::hclust(distMatrix, method = "average")
  
  # 2. Scale Clusters with Data Size (Forced Diversity)
  # Target roughly 1 domain per 3 terms, clamped between 5 and 25
  # This scales up complexity as the dataset gets richer.
  target_k <- max(5, min(25, ceiling(nTerms / 3)))
  target_k <- min(target_k, nTerms) # Don't exceed nTerms
  
  if (verbose) {
    message(sprintf("    -> Generating diverse palette: %d functional domains (Average Linkage).", target_k))
  }
  
  termDomains <- stats::cutree(hc, k = target_k)
  names(termDomains) <- terms
  return(list(termDomains = termDomains, numDomains = target_k))
}

#' @keywords internal
.assignDistinctColors <- function(numDomains) {
  if (numDomains <= 8) {
    palette <- RColorBrewer::brewer.pal(n = max(3, numDomains), name = "Dark2")[1:numDomains]
  } else if (numDomains <= 12) {
    palette <- RColorBrewer::brewer.pal(12, name = "Paired")[1:numDomains]
  } else {
    # Generate a larger diverse palette using colorspace if available
    if (requireNamespace("colorspace", quietly = TRUE)) {
      palette <- colorspace::qualitative_hcl(numDomains, palette = "Dark 3")
    } else {
      base <- RColorBrewer::brewer.pal(12, "Set3")
      palette <- grDevices::colorRampPalette(base)(numDomains)
    }
  }
  tibble::tibble(baseColor = palette)
}

#' @keywords internal
.blendColorsLAB <- function(df, minScore = 0) {
  valid <- !is.na(df$baseColor) & is.finite(df$score) & df$score > minScore
  
  if (!any(valid)) return("#CCCCCC")
  
  # Use LAB blending for better perceptual results if available
  if (requireNamespace("colorspace", quietly = TRUE)) {
    tryCatch({
      rgb_cols <- colorspace::hex2RGB(df$baseColor[valid])
      lab_cols <- as(rgb_cols, "LAB")
      weights <- df$score[valid] / sum(df$score[valid])
      # Weighted average in LAB space
      avg_L <- sum(lab_cols@coords[,1] * weights)
      avg_A <- sum(lab_cols@coords[,2] * weights)
      avg_B <- sum(lab_cols@coords[,3] * weights)
      return(colorspace::hex(colorspace::LAB(avg_L, avg_A, avg_B)))
    }, error = function(e) return("#CCCCCC"))
  } else {
    # Fallback to RGB mixing
    rgbMat <- t(grDevices::col2rgb(df$baseColor[valid]))
    weights <- df$score[valid] / sum(df$score[valid])
    res <- colSums(rgbMat * weights)
    return(grDevices::rgb(res[1], res[2], res[3], maxColorValue = 255))
  }
}