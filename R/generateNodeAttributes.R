utils::globalVariables(c("domainId", "Description", "p.adjust", "complexId",
                         "primaryFunctionalDomain", "topEnrichedFunctions",
                         "colorHex", "proteinCount"))

#' Generate Node Attributes for a Complex Network
#'
#' @description
#' Creates a detailed attribute table for each complex, suitable for network
#' visualization. Attributes include protein count, a primary functional
#' domain, top enriched functions, and a unique color representing the
#' complex's functional profile.
#'
#' @details
#' This function performs several steps to generate rich node attributes:
#' 1.  It aggregates all enriched terms from the input `enrichments` list.
#' 2.  A term-complex matrix is built, and terms are clustered based on their
#'     co-occurrence in complexes OR gene set overlap (if provided).
#' 3.  A unique color is assigned to each functional domain using a qualitative
#'     palette.
#' 4.  For each complex, it determines the primary functional domain based on
#'     the most significant enriched term.
#' 5.  A unique "blended" color is calculated for each complex by mixing the
#'     colors of its associated domains, weighted by significance.
#' 6.  Basic attributes like protein count and a list of proteins are also included.
#'
#' Complexes with no significant enrichments are assigned a default "Unenriched"
#' domain and a grey color.
#'
#' @param complexes A named list of protein complexes.
#' @param enrichments A named list of enrichment results, typically from
#'   `runComplexEnrichment`.
#' @param geneSetDb Optional named list where names are term IDs and values are
#'   character vectors of genes. If provided, term similarity will be calculated
#'   based on gene set overlap instead of co-occurrence. This results in more
#'   biologically meaningful color groupings.
#' @param similarityMethod The distance/similarity method. For gene set overlap:
#'   "overlap" (default), "jaccard", or "cosine". For co-occurrence: any method
#'   supported by `philentropy::distance`. Defaults to "overlap".
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return
#' A `tibble` where each row corresponds to a complex, with detailed attributes
#' for visualization.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette col2rgb rgb
#' @importFrom stats as.dist hclust cutree median setNames
#' @importFrom utils head
#' @importFrom methods as
#' @importFrom colorspace qualitative_hcl hex2RGB coords hex LAB
#'
#' @export
#' @examples
#' # --- Sample Data ---
#' complexes <- list(
#'   Cpx1 = c("A", "B", "C"),
#'   Cpx2 = c("C", "D", "E"),
#'   Cpx3 = c("F", "G")
#' )
#' enrichments <- list(
#'   Cpx1 = data.frame(ID = "GO:1", Description = "Term A", p.adjust = 0.01),
#'   Cpx2 = data.frame(ID = "GO:2", Description = "Term B", p.adjust = 0.02)
#' )
#'
#' # --- Without gene sets (co-occurrence) ---
#' nodeAttrs <- generateNodeAttributes(complexes, enrichments)
#'
#' # --- With gene sets (better functional grouping) ---
#' geneSets <- list(
#'   "GO:1" = c("A", "B", "X"),
#'   "GO:2" = c("A", "B", "Y"),  # High overlap with GO:1 -> similar color
#'   "GO:3" = c("M", "N", "O")   # No overlap -> different color
#' )
#' nodeAttrs <- generateNodeAttributes(complexes, enrichments,
#'                                     geneSetDb = geneSets)
#'
generateNodeAttributes <- function(complexes, enrichments,
                                   geneSetDb = NULL,
                                   similarityMethod = "overlap",
                                   verbose = TRUE) {
  if (verbose) {
    message("Generating core node attributes (function and color)...")
  }
  
  allTermsDf <- dplyr::bind_rows(enrichments, .id = "complexId")
  
  if (nrow(allTermsDf) == 0) {
    if (verbose) message("No enriched terms found; returning basic attributes.")
    proteinStr <-
      vapply(complexes, paste, collapse = ",", FUN.VALUE = character(1))
    return(
      tibble::tibble(
        complexId = names(complexes),
        proteinCount = lengths(complexes),
        proteins = proteinStr,
        primaryFunctionalDomain = "Unenriched",
        topEnrichedFunctions = NA_character_,
        colorHex = "#CCCCCC",
        sizeMapping = log2(pmax(1, lengths(complexes)))
      )
    )
  }
  
  enrichedComplexIds <- names(enrichments)
  allTerms <- unique(allTermsDf$ID)
  
  if (!is.null(geneSetDb)) {
    if (verbose) {
      message(sprintf(
        "    -> Clustering %d terms using gene set overlap (%s)",
        length(allTerms), similarityMethod
      ))
    }
    distMatrix <- .computeGeneSetDistance(allTerms, geneSetDb, similarityMethod)
  } else {
    if (verbose) {
      message(sprintf(
        "    -> Clustering %d terms using co-occurrence (%s)",
        length(allTerms), similarityMethod
      ))
    }
    distMatrix <- .computeCooccurrenceDistance(
      allTerms, enrichments, enrichedComplexIds, similarityMethod
    )
  }
  
  clusterResult <- .clusterTermsOptimal(distMatrix, allTerms, verbose)
  termDomains <- clusterResult$termDomains
  numDomains <- clusterResult$numDomains
  
  if (verbose) {
    message(sprintf("    -> Identified %d functional domains", numDomains))
  }
  
  termToDomainMap <-
    tibble::tibble(term = names(termDomains), domainId = termDomains)
  
  domainLabels <- allTermsDf %>%
    dplyr::left_join(termToDomainMap, by = c("ID" = "term")) %>%
    dplyr::group_by(domainId) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::summarise(
      domainLabel = dplyr::first(Description),
      .groups = "drop"
    )
  
  domainColors <- .assignDistinctColors(numDomains) %>%
    dplyr::mutate(domainId = seq_len(numDomains)) %>%
    dplyr::left_join(domainLabels, by = "domainId")
  
  allEnrichDf <- dplyr::bind_rows(
    lapply(names(enrichments), function(cid) {
      df <- enrichments[[cid]]
      if (!"p.adjust" %in% names(df)) df$p.adjust <- NA_real_
      df %>% dplyr::mutate(
        complexId = cid,
        score = -log10(pmax(p.adjust, 1e-300))
      )
    })
  ) %>%
    dplyr::left_join(termToDomainMap, by = c("ID" = "term")) %>%
    dplyr::left_join(domainColors, by = "domainId")
  
  complexSummary <- allEnrichDf %>%
    dplyr::group_by(complexId) %>%
    dplyr::group_map(.keep = TRUE, .f = function(df, key) {
      primary <- if (nrow(df) > 0 && any(!is.na(df$score))) {
        df$domainLabel[which.max(df$score)][1]
      } else NA_character_
      
      topFuncs <- if (nrow(df) > 0 && "Description" %in% names(df)) {
        paste(utils::head(df$Description[order(df$p.adjust)], 3),
              collapse = "; ")
      } else NA_character_
      
      colorHex <- .blendColorsLAB(df)
      
      tibble::tibble(
        complexId = unique(df$complexId),
        primaryFunctionalDomain = as.character(primary),
        topEnrichedFunctions = as.character(topFuncs),
        colorHex = colorHex
      )
    }) %>% dplyr::bind_rows()
  
  complexMeta <- tibble::tibble(
    complexId = names(complexes),
    proteinCount = lengths(complexes),
    proteins = vapply(complexes, paste, collapse = ",", FUN.VALUE = character(1))
  )
  
  finalDf <- complexMeta %>%
    dplyr::left_join(complexSummary, by = "complexId") %>%
    dplyr::mutate(
      primaryFunctionalDomain = dplyr::coalesce(
        primaryFunctionalDomain, "Unenriched"
      ),
      topEnrichedFunctions = dplyr::coalesce(
        topEnrichedFunctions, NA_character_
      ),
      colorHex = dplyr::coalesce(colorHex, "#CCCCCC"),
      sizeMapping = log2(pmax(1, proteinCount))
    )
  
  return(finalDf)
}


# === INTERNAL HELPER FUNCTIONS ===

#' @keywords internal
.computeGeneSetDistance <- function(terms, geneSetDb, method) {
  validTerms <- intersect(terms, names(geneSetDb))
  if (length(validTerms) == 0) {
    stop("No terms found in gene set database", call. = FALSE)
  }
  if (length(validTerms) < length(terms)) {
    warning(sprintf("%d terms not found in gene set database",
                    length(terms) - length(validTerms)))
  }
  termGeneSets <- geneSetDb[validTerms]
  nTerms <- length(termGeneSets)
  if (nTerms == 1) {
    return(stats::as.dist(matrix(0, 1, 1, dimnames = list(validTerms, validTerms))))
  }
  allGenes <- unique(unlist(termGeneSets, use.names = FALSE))
  geneTermMatrix <- matrix(0L, nrow = length(allGenes), ncol = nTerms,
                           dimnames = list(allGenes, validTerms))
  for (i in seq_along(termGeneSets)) {
    genes <- termGeneSets[[i]]
    geneTermMatrix[genes, i] <- 1L
  }
  intersectionMatrix <- t(geneTermMatrix) %*% geneTermMatrix
  setSizes <- colSums(geneTermMatrix)
  simMatrix <- matrix(0, nrow = nTerms, ncol = nTerms,
                      dimnames = list(validTerms, validTerms))
  diag(simMatrix) <- 1
  for (i in 1:(nTerms - 1)) {
    for (j in (i + 1):nTerms) {
      overlap <- intersectionMatrix[i, j]
      sim <- switch(method,
                    "overlap" = overlap / min(setSizes[i], setSizes[j]),
                    "jaccard" = overlap / (setSizes[i] + setSizes[j] - overlap),
                    "cosine" = overlap / sqrt(setSizes[i] * setSizes[j]),
                    overlap / min(setSizes[i], setSizes[j])
      )
      simMatrix[i, j] <- sim
      simMatrix[j, i] <- sim
    }
  }
  stats::as.dist(1 - simMatrix)
}

#' @keywords internal
.computeCooccurrenceDistance <- function(terms, enrichments, complexIds, method) {
  termComplexMatrix <- matrix(
    0L, nrow = length(terms), ncol = length(complexIds),
    dimnames = list(terms, complexIds)
  )
  idsList <- lapply(enrichments, function(x) x$ID)
  termIdx <- match(unlist(idsList, use.names = FALSE), terms)
  countsPerComplex <- lengths(idsList)
  complexIdx <- rep(seq_along(complexIds), times = countsPerComplex)
  if (length(termIdx) > 0) {
    termComplexMatrix[cbind(termIdx, complexIdx)] <- 1L
  }
  if (nrow(termComplexMatrix) == 1) {
    return(stats::as.dist(matrix(0, 1, 1, dimnames = list(terms, terms))))
  }
  if (method %in% c("jaccard", "overlap", "cosine")) {
    nTerms <- nrow(termComplexMatrix)
    simMatrix <- matrix(0, nrow = nTerms, ncol = nTerms,
                        dimnames = list(terms, terms))
    diag(simMatrix) <- 1
    for (i in 1:(nTerms - 1)) {
      for (j in (i + 1):nTerms) {
        vec_i <- termComplexMatrix[i, ]
        vec_j <- termComplexMatrix[j, ]
        overlap <- sum(vec_i & vec_j)
        sim <- switch(method,
                      "jaccard" = overlap / sum(vec_i | vec_j),
                      "cosine" = overlap / sqrt(sum(vec_i) * sum(vec_j)),
                      "overlap" = overlap / min(sum(vec_i), sum(vec_j)),
                      0
        )
        simMatrix[i, j] <- sim
        simMatrix[j, i] <- sim
      }
    }
    distMatrix <- stats::as.dist(1 - simMatrix)
  } else {
    if (!requireNamespace("philentropy", quietly = TRUE)) {
      stop("Package 'philentropy' required for method: ", method, call. = FALSE)
    }
    distMatrix <- philentropy::distance(
      termComplexMatrix, method = method, use.row.names = TRUE
    )
  }
  return(distMatrix)
}

#' @keywords internal
.clusterTermsOptimal <- function(distMatrix, terms, verbose) {
  nTerms <- length(terms)
  if (nTerms <= 1) {
    return(list(termDomains = stats::setNames(1, terms), numDomains = 1))
  }
  hc <- stats::hclust(distMatrix, method = "ward.D2")
  maxK <- min(15, nTerms - 1, 20)
  minK <- 2
  if (maxK > minK && nTerms > 3) {
    wss <- sapply(minK:min(maxK, 10), function(k) {
      clusters <- stats::cutree(hc, k = k)
      sum(tapply(seq_along(clusters), clusters, function(idx) {
        if (length(idx) <= 1) return(0)
        subDist <- as.matrix(distMatrix)[idx, idx]
        sum(subDist^2) / (2 * length(idx))
      }))
    })
    if (length(wss) > 2) {
      drops <- -diff(wss)
      numDomains <- (minK:min(maxK, 10))[which.max(drops) + 1]
    } else {
      numDomains <- minK
    }
    numDomains <- min(numDomains, max(8, ceiling(nTerms / 5)))
  } else {
    numDomains <- min(minK, nTerms)
  }
  termDomains <- stats::cutree(hc, k = numDomains)
  names(termDomains) <- terms
  return(list(termDomains = termDomains, numDomains = numDomains))
}

#' @keywords internal
.assignDistinctColors <- function(numDomains) {
  if (numDomains <= 8) {
    palette <- RColorBrewer::brewer.pal(n = max(3, numDomains), name = "Dark2")
    if (numDomains < 3) palette <- palette[1:numDomains]
  } else if (numDomains <= 12) {
    palette <- RColorBrewer::brewer.pal(12, name = "Paired")
  } else {
    if (requireNamespace("colorspace", quietly = TRUE)) {
      palette <- colorspace::qualitative_hcl(numDomains, palette = "Dark 3")
    } else {
      basePalette <- RColorBrewer::brewer.pal(9, "Set1")
      palette <- grDevices::colorRampPalette(basePalette)(numDomains)
    }
  }
  tibble::tibble(baseColor = palette)
}

#' @keywords internal
.blendColorsLAB <- function(df, minScore = -log10(0.05)) {
  valid <- !is.na(df$baseColor) &
    !is.na(df$score) &
    is.finite(df$score) &
    df$score > minScore
  if (!any(valid) || sum(df$score[valid], na.rm = TRUE) == 0) {
    return("#CCCCCC")
  }
  if (requireNamespace("colorspace", quietly = TRUE)) {
    tryCatch({
      # CORRECTED LINE: Use the public S4 coercion method 'as(..., "LAB")'
      labColors <- t(sapply(df$baseColor[valid], function(col) {
        colorspace::coords(methods::as(colorspace::hex2RGB(col), "LAB"))
      }))
      weights <- df$score[valid] / sum(df$score[valid])
      labMean <- colSums(labColors * weights)
      return(colorspace::hex(colorspace::LAB(labMean[1], labMean[2], labMean[3])))
    }, error = function(e) {
      # Fall through
    })
  }
  rgbMat <- t(grDevices::col2rgb(df$baseColor[valid, drop = TRUE]))
  weights <- df$score[valid] / sum(df$score[valid])
  rgbVal <- colSums(rgbMat * weights)
  grDevices::rgb(rgbVal[1], rgbVal[2], rgbVal[3], maxColorValue = 255)
}