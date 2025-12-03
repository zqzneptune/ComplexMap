utils::globalVariables(c("Term", "Count", "Ratio", "Fold", "PopHits", "TermSize"))

#' Run Gene Set Enrichment Analysis on a List of Complexes
#'
#' @description
#' Performs enrichment analysis for each complex against a GMT background.
#'
#' @details
#' **Systems Biology Rationale:**
#' To support the generation of a specific and diverse landscape, this function
#' calculates and returns **Fold Enrichment** and **Term Size** in addition to
#' standard p-values.
#' 
#' - **Fold Enrichment** is used by downstream functions to prioritize specific
#'   biological labels over generic ones.
#' - **Term Size** (PopHits) allows you to filter out overly broad terms
#'   (e.g., those containing > 10% of the genome) if desired.
#'
#' @param complexList A named list of protein complexes (character vectors).
#' @param gmt A named list of gene sets (character vectors).
#' @param pAdjustMethod Method for multiple testing correction ("Benjamini", "Bonferroni", "FDR").
#'   Defaults to "Benjamini".
#' @param pValueCutoff Adjusted p-value cutoff. Defaults to 0.05.
#' @param verbose Logical.
#'
#' @return
#' A named list of tibbles. Each tibble contains:
#' - `ID`: Term ID/Name
#' - `Description`: Term Description
#' - `p.adjust`: Adjusted p-value
#' - `Fold`: Fold Enrichment (Observed/Expected)
#' - `TermSize`: Number of genes in the term (Background)
#' - `Count`: Number of genes in the term (Overlap)
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
runComplexEnrichment <- function(complexList,
                                 gmt,
                                 pAdjustMethod="Benjamini",
                                 pValueCutoff=0.05,
                                 verbose=TRUE) {
  if (verbose) {
    message(sprintf("Running enrichment for %d complexes (Cutoff: %.2f)...", 
                    length(complexList), pValueCutoff))
  }
  
  enrichmentResults <- list()
  
  for (complexId in names(complexList)) {
    complexGenes <- complexList[[complexId]]
    result <- .runHypergeometricTest(geneSet = complexGenes, gmt = gmt)
    
    if (!is.null(result) && nrow(result) > 0) {
      validMethods <- c("Benjamini", "Bonferroni", "FDR")
      if (!pAdjustMethod %in% validMethods) {
        stop(sprintf("pAdjustMethod '%s' not found.", pAdjustMethod))
      }
      
      # Filter by significance
      significantResult <- result %>%
        dplyr::filter(!!rlang::sym(pAdjustMethod) < pValueCutoff)
      
      if (nrow(significantResult) > 0) {
        # Select and Rename columns for the standardized ComplexMap format
        enrichmentResults[[complexId]] <- significantResult %>%
          dplyr::select(
            ID = Term,
            Description = Term,
            p.adjust = dplyr::all_of(pAdjustMethod),
            Fold,
            TermSize = PopHits, # Renamed for clarity: Size of pathway in Universe
            Count,
            Ratio
          ) %>%
          dplyr::arrange(p.adjust) # Sort by significance by default
      }
    }
  }
  
  if (verbose) {
    message(sprintf("Annotation complete. Found significant terms for %d complexes.",
                    length(enrichmentResults)))
  }
  return(enrichmentResults)
}

#' @keywords internal
.runHypergeometricTest <- function(geneSet, gmt) {
  # 1. Define Universe (Background)
  gsUniverse <- unique.default(unlist(gmt, use.names = FALSE))
  popTotal <- length(gsUniverse)
  
  # 2. Intersect Input with Universe
  # Genes in the complex that are NOT in the GMT/Universe are excluded from the test
  geneSetInUniverse <- intersect(geneSet, gsUniverse)
  listTotal <- length(geneSetInUniverse)
  
  if (listTotal == 0) return(NULL)
  
  # 3. Identify overlaps
  # Logic: Which terms in GMT overlap with our filtered gene set?
  hasOverlap <- vapply(gmt, function(termGenes) {
    # Any shared genes?
    length(intersect(termGenes, geneSetInUniverse)) > 0
  }, logical(1))
  
  gmtFiltered <- gmt[hasOverlap]
  if (length(gmtFiltered) == 0) return(NULL)
  
  # 4. Calculate Statistics
  popHits <- lengths(gmtFiltered)
  count <- lengths(lapply(gmtFiltered, intersect, geneSetInUniverse))
  
  # Expected Count = (ListSize * TermSize) / UniverseSize
  expected <- (listTotal * popHits) / popTotal
  
  # Fold Enrichment = Observed / Expected
  fold <- count / expected
  
  # Ratio = Count / ListSize (How much of the complex is covered?)
  ratio <- count / listTotal
  
  # 5. Hypergeometric Test (Phyper)
  # q = count - 1 (because lower.tail = FALSE is P[X > x])
  # m = popHits (white balls)
  # n = popTotal - popHits (black balls)
  # k = listTotal (draws)
  pValue <- stats::phyper(
    count - 1, popHits, popTotal - popHits, listTotal, lower.tail = FALSE
  )
  
  data.frame(
    Term = names(gmtFiltered),
    Count = count,
    Ratio = ratio,
    PValue = pValue,
    Fold = fold,
    ListTotal = listTotal,
    PopHits = popHits,
    PopTotal = popTotal,
    Benjamini = stats::p.adjust(pValue, method = "BH"),
    Bonferroni = stats::p.adjust(pValue, method = "bonferroni"),
    FDR = stats::p.adjust(pValue, method = "fdr"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}