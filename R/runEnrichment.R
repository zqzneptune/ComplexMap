utils::globalVariables(c("Term", "Count", "Ratio", "Fold"))
#' Run Gene Set Enrichment Analysis on a List of Complexes
#'
#' @description
#' This function performs a hypergeometric-based enrichment analysis for each
#' protein complex in a list against a provided gene set matrix (GMT).
#'
#' @details
#' For each complex, this function calculates the over-representation of
#' functional terms (e.g., GO terms, pathways) from the GMT file. It uses a
#' hypergeometric test to compute a p-value, which is then adjusted for
#' multiple testing.
#'
#' Only the terms that are significant after filtering by the `pValueCutoff`
#' are retained in the final output.
#'
#' @param complexList A named list where each element is a character vector of
#'   gene/protein identifiers representing a complex.
#' @param gmt A named list where each element is a character vector of genes,
#'   representing a functional gene set (e.g., from a GMT file).
#' @param pAdjustMethod A character string specifying the p-value adjustment
#'   method to use for filtering. Must be one of "Benjamini", "Bonferroni",
#'   or "FDR". Defaults to "Benjamini".
#' @param pValueCutoff A numeric value used as the cutoff for significance on
#'   the adjusted p-value. Defaults to 0.05.
#' @param verbose A logical value indicating whether to print progress messages.
#'   Defaults to `TRUE`.
#'
#' @return
#' A named list where each name corresponds to a `complexId` from the input.
#' Each element is a data frame containing the significant enrichment results
#' for that complex, with columns: `ID`, `Description`, `p.adjust`, `Count`,
#' `Ratio`, and `Fold`.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @examples
#' # --- Sample Data ---
#' # 1. Complexes to be tested
#' complex1 <- c("A", "B", "C", "D")
#' complex2 <- c("F", "G", "H")
#' myComplexes <- list(Cpx1 = complex1, Cpx2 = complex2)
#'
#' # 2. Gene Set Matrix (e.g., GO terms or pathways)
#' term1 <- c("A", "B", "C", "X", "Y") # Enriched in Cpx1
#' term2 <- c("F", "G", "Z")          # Enriched in Cpx2
#' term3 <- c("L", "M", "N")          # Not enriched
#' myGmt <- list(Term1 = term1, Term2 = term2, Term3 = term3)
#'
#' # --- Run Enrichment ---
#' enrichment <- runComplexEnrichment(myComplexes, myGmt)
#' print(enrichment)
#'
#' @export
#' @importFrom dplyr filter all_of
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @importFrom stats phyper p.adjust
#'
runComplexEnrichment <- function(complexList,
                                 gmt,
                                 pAdjustMethod = "Benjamini",
                                 pValueCutoff = 0.05,
                                 verbose = TRUE) {
  
  if (verbose) {
    message(
      sprintf("Running enrichment for %d complexes...", length(complexList))
    )
  }
  
  enrichmentResults <- list()
  
  for (complexId in names(complexList)) {
    complexGenes <- complexList[[complexId]]
    
    # Use the internal helper function for the core calculation
    result <- .runHypergeometricTest(geneSet = complexGenes, gmt = gmt)
    
    if (!is.null(result) && nrow(result) > 0) {
      # Check if the chosen p-adjust method exists as a column
      validMethods <- c("Benjamini", "Bonferroni", "FDR")
      if (!pAdjustMethod %in% validMethods) {
        stop(
          sprintf("pAdjustMethod '%s' not found. Choose one of: %s",
                  pAdjustMethod, paste(validMethods, collapse = ", "))
        )
      }
      
      # Filter for significant results using dplyr and NSE
      significantResult <- result %>%
        filter(!!sym(pAdjustMethod) < pValueCutoff)
      
      if (nrow(significantResult) > 0) {
        # Select and rename columns for a clean output format
        enrichmentResults[[complexId]] <- significantResult %>%
          dplyr::select(
            ID = Term,
            Description = Term,
            p.adjust = all_of(pAdjustMethod),
            Count,
            Ratio,
            Fold
          )
      }
    }
  }
  
  if (verbose) {
    message(
      sprintf("Annotation complete. Found terms for %d complexes.",
              length(enrichmentResults))
    )
  }
  return(enrichmentResults)
}

# Internal helper to perform a hypergeometric test for one gene set.
# Not exported.
#
# @param geneSet A character vector of genes in the query set (e.g., complex).
# @param gmt A named list representing the gene set universe.
# @return A data.frame with raw enrichment statistics, or NULL.
.runHypergeometricTest <- function(geneSet, gmt) {
  
  gsUniverse <- unique.default(unlist(gmt, use.names = FALSE))
  popTotal <- length(gsUniverse)
  listTotal <- length(intersect(geneSet, gsUniverse))
  
  if (listTotal == 0) {
    return(NULL) # No overlap between the complex and the universe
  }
  
  # Filter GMT to include only terms with at least one overlapping gene
  hasOverlap <- vapply(gmt,
                       function(termGenes) length(intersect(termGenes, geneSet)) > 0,
                       logical(1)
  )
  gmtFiltered <- gmt[hasOverlap]
  
  if (length(gmtFiltered) == 0) {
    return(NULL) # No gene sets found with overlapping genes
  }
  
  # Calculate enrichment statistics using vectorized operations
  popHits <- lengths(gmtFiltered) # Size of each term in the universe
  count <- lengths(lapply(gmtFiltered, intersect, geneSet)) # Overlap size
  ratio <- count / listTotal
  fold <- (count / listTotal) / (popHits / popTotal)
  
  # Hypergeometric test p-value: P(X >= count)
  pValue <- phyper(
    count - 1, popHits, popTotal - popHits, listTotal, lower.tail = FALSE
  )
  
  # Pre-calculate all common adjusted p-values
  df <- data.frame(
    Term = names(gmtFiltered),
    Count = count,
    Ratio = ratio,
    PValue = pValue,
    Fold = fold,
    ListTotal = listTotal,
    PopHits = popHits,
    PopTotal = popTotal,
    Benjamini = p.adjust(pValue, method = "BH"),
    Bonferroni = p.adjust(pValue, method = "bonferroni"),
    FDR = p.adjust(pValue, method = "fdr"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(df)
}