# Run Gene Set Enrichment Analysis on a List of Complexes

Performs enrichment analysis for each complex against a GMT background.

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

  A named list of protein complexes (character vectors).

- gmt:

  A named list of gene sets (character vectors).

- pAdjustMethod:

  Method for multiple testing correction ("Benjamini", "Bonferroni",
  "FDR"). Defaults to "Benjamini".

- pValueCutoff:

  Adjusted p-value cutoff. Defaults to 0.05.

- verbose:

  Logical.

## Value

A named list of tibbles. Each tibble contains: - \`ID\`: Term ID/Name -
\`Description\`: Term Description - \`p.adjust\`: Adjusted p-value -
\`Fold\`: Fold Enrichment (Observed/Expected) - \`TermSize\`: Number of
genes in the term (Background) - \`Count\`: Number of genes in the term
(Overlap)

## Details

\*\*Systems Biology Rationale:\*\* To support the generation of a
specific and diverse landscape, this function calculates and returns
\*\*Fold Enrichment\*\* and \*\*Term Size\*\* in addition to standard
p-values.

\- \*\*Fold Enrichment\*\* is used by downstream functions to prioritize
specific biological labels over generic ones. - \*\*Term Size\*\*
(PopHits) allows you to filter out overly broad terms (e.g., those
containing \> 10

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
