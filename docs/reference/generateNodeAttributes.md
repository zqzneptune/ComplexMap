# Generate Node Attributes for a Complex Network

Creates a detailed attribute table for each complex, emphasizing
functional specificity for network visualization, and assigns
professional HCL colors.

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

  A named list of enrichment results. Must contain 'Fold'.

- geneSetDb:

  Optional named list of gene sets for semantic clustering.

- similarityMethod:

  Distance method for clustering ("jaccard", etc.).

- verbose:

  Logical.

## Value

A \`tibble\` with complex attributes.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
