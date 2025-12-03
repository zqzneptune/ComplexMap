# Generate Node Attributes for a Complex Network

Creates a detailed attribute table for each complex, emphasizing
functional specificity for network visualization.

## Usage

``` r
generateNodeAttributes(
  complexes,
  enrichments,
  geneSetDb = NULL,
  similarityMethod = "jaccard",
  verbose = TRUE
)
```

## Arguments

- complexes:

  A named list of protein complexes.

- enrichments:

  A named list of enrichment results. Must contain 'Fold' column for
  specificity weighting.

- geneSetDb:

  Optional named list of gene sets for semantic clustering.

- similarityMethod:

  Distance method for clustering ("jaccard", "overlap", etc.). Defaults
  to "jaccard" to penalize size differences.

- verbose:

  Logical.

## Value

A \`tibble\` with complex attributes.

## Details

This function performs several steps to generate rich node
attributes: 1. It aggregates all enriched terms from the input
\`enrichments\` list. 2. It calculates a \*\*Specificity Score\*\* for
each term using the formula: \`-log10(p.adjust) \* log2(Fold)\`. This
ensures that specific, high-fold enrichment terms are prioritized over
broad, generic terms. 3. A term-complex matrix is built, and terms are
clustered. The clustering uses \*\*Average Linkage\*\* and forces a
higher number of clusters to preserve functional diversity (avoiding
"monochromatic" maps). 4. A unique color is assigned to each functional
domain. 5. For each complex, the \*\*Primary Functional Domain\*\* is
assigned to the enriched term with the highest Specificity Score. 6. A
unique "blended" color is calculated by mixing domain colors, weighted
by the Specificity Score.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
