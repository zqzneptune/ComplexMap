# Run Gene Set Enrichment Analysis on a List of Complexes

This function performs a hypergeometric-based enrichment analysis for
each protein complex in a list against a provided gene set matrix (GMT).

## Usage

``` r
runComplexEnrichment(
  complexList,
  gmt,
  pAdjustMethod = "Benjamini",
  pValueCutoff = 0.05,
  verbose = TRUE
)
```

## Arguments

- complexList:

  A named list where each element is a character vector of gene/protein
  identifiers representing a complex.

- gmt:

  A named list where each element is a character vector of genes,
  representing a functional gene set (e.g., from a GMT file).

- pAdjustMethod:

  A character string specifying the p-value adjustment method to use for
  filtering. Must be one of "Benjamini", "Bonferroni", or "FDR".
  Defaults to "Benjamini".

- pValueCutoff:

  A numeric value used as the cutoff for significance on the adjusted
  p-value. Defaults to 0.05.

- verbose:

  A logical value indicating whether to print progress messages.
  Defaults to \`TRUE\`.

## Value

A named list where each name corresponds to a \`complexId\` from the
input. Each element is a data frame containing the significant
enrichment results for that complex, with columns: \`ID\`,
\`Description\`, \`p.adjust\`, \`Count\`, \`Ratio\`, and \`Fold\`.

## Details

For each complex, this function calculates the over-representation of
functional terms (e.g., GO terms, pathways) from the GMT file. It uses a
hypergeometric test to compute a p-value, which is then adjusted for
multiple testing.

Only the terms that are significant after filtering by the
\`pValueCutoff\` are retained in the final output.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data ---
# 1. Complexes to be tested
complex1 <- c("A", "B", "C", "D")
complex2 <- c("F", "G", "H")
myComplexes <- list(Cpx1 = complex1, Cpx2 = complex2)

# 2. Gene Set Matrix (e.g., GO terms or pathways)
term1 <- c("A", "B", "C", "X", "Y") # Enriched in Cpx1
term2 <- c("F", "G", "Z")          # Enriched in Cpx2
term3 <- c("L", "M", "N")          # Not enriched
myGmt <- list(Term1 = term1, Term2 = term2, Term3 = term3)

# --- Run Enrichment ---
enrichment <- runComplexEnrichment(myComplexes, myGmt)
#> Running enrichment for 2 complexes...
#> Annotation complete. Found terms for 0 complexes.
print(enrichment)
#> list()
```
