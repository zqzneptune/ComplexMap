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
#'     co-occurrence in complexes using the specified similarity metric. This
#'     groups related functional terms into "functional domains".
#' 3.  A unique color is assigned to each functional domain using a qualitative
#'     palette from `RColorBrewer`.
#' 4.  For each complex, it determines the primary functional domain based on
#'     the most significant enriched term.
#' 5.  A unique "blended" color is calculated for each complex by mixing the
#'     colors of its associated domains, weighted by the significance
#'     (-log10 p-value) of the enriched terms.
#' 6.  Basic attributes like protein count and a list of proteins are also included.
#'
#' Complexes with no significant enrichments are assigned a default "Unenriched"
#' domain and a grey color.
#'
#' @param complexes A named list of protein complexes.
#' @param enrichments A named list of enrichment results, typically from
#'   `runComplexEnrichment`. Each element should be a data frame with at least
#'   `ID`, `Description`, and `p.adjust` columns.
#' @param similarityMethod The distance/similarity method passed to
#'   `philentropy::distance` for clustering terms. Defaults to "jaccard".
#' @param verbose A logical value indicating whether to print progress messages.
#'   Defaults to `TRUE`.
#'
#' @return
#' A `tibble` where each row corresponds to a complex. The columns include:
#' `complexId`, `proteinCount`, `proteins`, `primaryFunctionalDomain`,
#' `topEnrichedFunctions`, `colorHex`, and `sizeMapping`.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @examples
#' # --- Sample Data ---
#' complexes <- list(
#'   Cpx1 = c("A", "B", "C"),
#'   Cpx2 = c("C", "D", "E"),
#'   Cpx3 = c("F", "G") # Unenriched
#' )
#' enrichments <- list(
#'   Cpx1 = data.frame(ID = "GO:1", Description = "Term A", p.adjust = 0.01),
#'   Cpx2 = data.frame(ID = "GO:2", Description = "Term B", p.adjust = 0.02)
#' )
#'
#' # --- Generate Node Attributes ---
#' nodeAttrs <- generateNodeAttributes(complexes, enrichments)
#' print(nodeAttrs)
#'
#' @export
#'
generateNodeAttributes <- function(complexes, enrichments,
                                   similarityMethod="jaccard",
                                   verbose=TRUE) {
  if (verbose) {
    message("Generating core node attributes (function and color)...")
    message(
      sprintf(
        "    -> Clustering terms using '%s' similarity.",
        similarityMethod
      )
    )
  }
  
  allTermsDf <- dplyr::bind_rows(enrichments, .id="complexId")
  
  # Handle case with no enriched terms across all complexes
  if (nrow(allTermsDf) == 0) {
    if (verbose) message("No enriched terms found; returning basic attributes.")
    proteinStr <- 
      vapply(complexes, paste, collapse=",", FUN.VALUE=character(1))
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
  
  # --- 1. Cluster functional terms into domains ---
  enrichedComplexIds <- names(enrichments)
  allTerms <- unique(allTermsDf$ID)
  termComplexMatrix <- matrix(
    0L, nrow=length(allTerms), ncol=length(enrichedComplexIds),
    dimnames=list(allTerms, enrichedComplexIds)
  )
  
  # Efficiently populate the binary term-complex matrix
  idsList <- lapply(enrichments, function(x) x$ID)
  termIdx <- match(unlist(idsList, use.names=FALSE), allTerms)
  countsPerComplex <- lengths(idsList)
  complexIdx <- rep(seq_along(enrichedComplexIds), times=countsPerComplex)
  if (length(termIdx) > 0) termComplexMatrix[cbind(termIdx, complexIdx)] <- 1L
  
  if (nrow(termComplexMatrix) > 1) {
    if (!requireNamespace("philentropy", quietly=TRUE)) {
      stop("Package 'philentropy' is required for term clustering.",
           call.=FALSE)
    }
    distMatrix <- philentropy::distance(
      termComplexMatrix, method=similarityMethod, use.row.names=TRUE
    )
    
    # Check if there is more than one unique row to cluster
    if (nrow(distMatrix) > 1 && any(distMatrix[lower.tri(distMatrix)] > 0)) {
      hc <- stats::hclust(stats::as.dist(distMatrix), method="average")
      numDomains <- min(15, nrow(distMatrix))
      termDomains <- stats::cutree(hc, k=numDomains)
    } else {
      numDomains <- nrow(termComplexMatrix)
      termDomains <- stats::setNames(
        seq_len(numDomains), rownames(termComplexMatrix)
      )
    }
  } else {
    numDomains <- 1
    termDomains <- stats::setNames(1, rownames(termComplexMatrix))
  }
  
  termToDomainMap <- 
    tibble::tibble(term=names(termDomains), domainId=termDomains)
  
  # --- 2. Create domain labels and color palette ---
  domainLabels <- allTermsDf %>%
    dplyr::left_join(termToDomainMap, by=c("ID"="term")) %>%
    dplyr::group_by(domainId) %>%
    dplyr::summarise(
      domainLabel=if (any(!is.na(Description))) Description[1] else NA_character_,
      .groups="drop"
    )
  
  basePalette <- RColorBrewer::brewer.pal(n=9, name="Set1")
  if (numDomains > length(basePalette)) {
    palette <- grDevices::colorRampPalette(basePalette)(numDomains)
  } else {
    palette <- basePalette[seq_len(numDomains)]
  }
  domainColors <- tibble::tibble(
    domainId = seq_len(numDomains), baseColor = palette
  ) %>%
    dplyr::left_join(domainLabels, by="domainId")
  
  # --- 3. Calculate per-complex attributes ---
  allEnrichDf <- dplyr::bind_rows(
    lapply(names(enrichments), function(cid) {
      df <- enrichments[[cid]]
      if (!"p.adjust" %in% names(df)) df$p.adjust <- NA_real_
      df %>% dplyr::mutate(complexId=cid, score=-log10(p.adjust))
    })
  ) %>%
    dplyr::left_join(termToDomainMap, by=c("ID"="term")) %>%
    dplyr::left_join(domainColors, by="domainId")
  
  complexSummary <- allEnrichDf %>%
    dplyr::group_by(complexId) %>%
    dplyr::group_map(.keep=TRUE, .f=function(df, key) {
      primary <- if (nrow(df) > 0 && any(!is.na(df$score))) {
        df$domainLabel[which.max(df$score)][1]
      } else NA_character_
      topFuncs <- if (nrow(df) > 0 && "Description" %in% names(df)) {
        paste(utils::head(df$Description[order(df$p.adjust)], 3),
              collapse=", ")
      } else NA_character_
      
      colorHex <- "#CCCCCC"
      valid <- !is.na(df$baseColor) & !is.na(df$score) & is.finite(df$score)
      if (any(valid) && sum(df$score[valid], na.rm=TRUE) > 0) {
        rgbMat <- t(grDevices::col2rgb(df$baseColor[valid, drop=TRUE]))
        weights <- df$score[valid] / sum(df$score[valid], na.rm=TRUE)
        rgbVal <- colSums(rgbMat * weights, na.rm=TRUE)
        colorHex <- grDevices::rgb(
          rgbVal[1], rgbVal[2], rgbVal[3], maxColorValue=255
        )
      }
      tibble::tibble(
        complexId = unique(df$complexId),
        primaryFunctionalDomain = as.character(primary),
        topEnrichedFunctions = as.character(topFuncs),
        colorHex = colorHex
      )
    }) %>% dplyr::bind_rows()
  
  # --- 4. Combine with base attributes and finalize ---
  complexMeta <- tibble::tibble(
    complexId = names(complexes),
    proteinCount = lengths(complexes),
    proteins = vapply(complexes, paste, collapse=",", FUN.VALUE=character(1))
  )
  
  finalDf <- complexMeta %>%
    dplyr::left_join(complexSummary, by="complexId") %>%
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