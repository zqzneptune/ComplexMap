utils::globalVariables(c("domainId", "Description", "p.adjust", "complexId",
                         "primaryFunctionalDomain", "topEnrichedFunctions",
                         "colorHex", "proteinCount", "Fold", "score", "as",
                         "domainLabel", "topIdx", "baseColor", "head"))

#' Generate Node Attributes for a Complex Network
#'
#' @description
#' Creates a detailed attribute table for each complex, emphasizing functional
#' specificity for network visualization, and assigns professional HCL colors.
#'
#' @param complexes A named list of protein complexes.
#' @param enrichments A named list of enrichment results. Must contain 'Fold'.
#' @param geneSetDb Optional named list of gene sets for semantic clustering.
#' @param similarityMethod Distance method for clustering ("jaccard", etc.).
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
    message("Generating node attributes (prioritizing functional specificity colors)...")
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
      colorHex = "#D3D3D3"
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
  
  clusterResult <- .clusterTermsOptimal(distMatrix, allTerms, verbose)
  termDomains <- clusterResult$termDomains
  numDomains <- clusterResult$numDomains
  
  termToDomainMap <- tibble::tibble(term = names(termDomains), domainId = termDomains)
  
  # --- 2. Build Semantic Palette ---
  # We compute the distance between CLUSTERS to order the colors
  # This ensures the "color wheel" matches the "functional wheel"
  domainDist <- .computeDomainDistance(distMatrix, termDomains)
  domainColors <- .assignProfessionalPalette(numDomains, domainDist)
  
  # --- 3. Create Domain Labels ---
  domainLabels <- allTermsDf %>%
    dplyr::left_join(termToDomainMap, by = c("ID" = "term")) %>%
    dplyr::group_by(domainId) %>%
    dplyr::summarise(
      domainLabel = Description[which.max(Fold)][1], 
      .groups = "drop"
    ) %>%
    dplyr::left_join(domainColors, by = "domainId")
  
  # --- 4. Calculate Per-Complex Attributes (Weighted by Specificity) ---
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
    dplyr::left_join(domainLabels, by = "domainId")
  
  complexSummary <- allEnrichDf %>%
    dplyr::group_by(complexId) %>%
    dplyr::summarise(
      topIdx = which.max(score)[1],
      primaryFunctionalDomain = domainLabel[topIdx],
      colorHex = baseColor[topIdx],
      topEnrichedFunctions = paste(unique(utils::head(Description[order(score, decreasing = TRUE)], 3)), collapse = "; "),
      .groups = "drop"
    )
  
  # --- 5. Finalize ---
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
      colorHex = dplyr::coalesce(colorHex, "#D3D3D3")
    )
  
  return(finalDf)
}


# ==============================================================================
# === INTERNAL HELPER FUNCTIONS ===
# ==============================================================================

#' @keywords internal
.computeDomainDistance <- function(termDist, termDomains) {
  if (is.null(termDist)) return(NULL)
  mat <- as.matrix(termDist)
  numDomains <- max(termDomains)
  domainDist <- matrix(0, nrow = numDomains, ncol = numDomains)
  
  for (i in seq_len(numDomains)) {
    for (j in seq_len(numDomains)) {
      if (i == j) {
        domainDist[i, j] <- 0
      } else {
        idx_i <- which(termDomains == i)
        idx_j <- which(termDomains == j)
        if (length(idx_i) > 0 && length(idx_j) > 0) {
          d <- mean(mat[idx_i, idx_j, drop = FALSE])
          domainDist[i, j] <- d
          domainDist[j, i] <- d
        }
      }
    }
  }
  return(stats::as.dist(domainDist))
}

#' Assign a Professional Color Palette using RColorBrewer
#'
#' @description
#' Generates a professional, semantically-ordered color palette.
#' Dynamically scales from `RColorBrewer` palettes (Dark2, Paired, Set3) based
#' on the number of domains, and maps functionally similar clusters to 
#' adjacent colors.
#'
#' @param numDomains Integer. The number of unique functional domains (clusters) to color.
#' @param distMatrix A `dist` object representing the distance between domains. Optional.
#'
#' @return A `tibble` containing two columns: `domainId` and `baseColor` (hex codes).
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats cmdscale
#' @importFrom tibble tibble
#'
#' @keywords internal
#' @noRd
.assignProfessionalPalette <- function(numDomains, distMatrix = NULL) {
  
  # 1. Generate base colors using RColorBrewer
  if (numDomains <= 8) {
    # Dark2 is highly professional and colorblind-friendly for small sets
    palette <- RColorBrewer::brewer.pal(n = max(3, numDomains), name = "Dark2")[seq_len(numDomains)]
  } else if (numDomains <= 12) {
    # Paired gives beautiful contrasts for medium-sized sets
    palette <- RColorBrewer::brewer.pal(12, name = "Paired")[seq_len(numDomains)]
  } else {
    # Set3 is pastel and visually distinct; extend it smoothly if > 12 domains
    base_palette <- RColorBrewer::brewer.pal(12, "Set3")
    palette <- grDevices::colorRampPalette(base_palette)(numDomains)
  }
  
  # 2. Semantic Ordering (Optional but recommended)
  # Ensure functionally similar domains are assigned adjacent colors
  color_order <- seq_len(numDomains)
  
  if (!is.null(distMatrix) && numDomains > 2) {
    try({
      # 1D Multidimensional Scaling maps similar clusters to nearby numeric values
      fit <- stats::cmdscale(distMatrix, k = 1)
      color_order <- order(fit)
    }, silent = TRUE)
  }
  
  # Apply the sorted order to the palette
  palette <- palette[color_order]
  
  return(tibble::tibble(
    domainId = seq_len(numDomains),
    baseColor = palette
  ))
}

#' @keywords internal
.computeGeneSetDistance <- function(terms, geneSetDb, method) {
  validTerms <- intersect(terms, names(geneSetDb))
  if (length(validTerms) == 0) stop("No terms found in gene set database")
  
  termGeneSets <- geneSetDb[validTerms]
  nTerms <- length(termGeneSets)
  if (nTerms == 1) return(stats::as.dist(matrix(0, 1, 1, dimnames = list(validTerms, validTerms))))
  
  allGenes <- unique(unlist(termGeneSets, use.names = FALSE))
  geneTermMatrix <- Matrix::sparseMatrix(
    i = match(unlist(termGeneSets), allGenes),
    j = rep(seq_along(termGeneSets), lengths(termGeneSets)),
    x = 1, dims = c(length(allGenes), nTerms),
    dimnames = list(allGenes, validTerms)
  )
  
  intersectionMatrix <- as.matrix(Matrix::crossprod(geneTermMatrix)) # T(M) %*% M
  setSizes <- Matrix::colSums(geneTermMatrix)
  
  simMatrix <- matrix(0, nrow = nTerms, ncol = nTerms, dimnames = list(validTerms, validTerms))
  
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
    raw_dist <- philentropy::distance(termComplexMatrix, method = method, use.row.names = TRUE)
    return(stats::as.dist(raw_dist))
  }
  
  return(stats::dist(termComplexMatrix, method = "binary"))
}

#' @keywords internal
.clusterTermsOptimal <- function(distMatrix, terms, verbose) {
  nTerms <- length(terms)
  if (nTerms <= 1) {
    return(list(termDomains = stats::setNames(1, terms), numDomains = 1))
  }
  
  hc <- stats::hclust(distMatrix, method = "average")
  
  target_k <- max(5, min(25, ceiling(nTerms / 3)))
  target_k <- min(target_k, nTerms) 
  
  if (verbose) {
    message(sprintf("    -> Generating diverse palette: %d functional domains (Average Linkage).", target_k))
  }
  
  termDomains <- stats::cutree(hc, k = target_k)
  names(termDomains) <- terms
  return(list(termDomains = termDomains, numDomains = target_k))
}