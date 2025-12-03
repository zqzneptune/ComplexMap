# Build a Complex-Complex Interaction Network

Constructs a network of complexes where edges represent similarity.

## Usage

``` r
buildComplexNetwork(
  complexes,
  enrichments,
  mode = "combined",
  similarityMethod = "jaccard",
  alpha = 0.75,
  nCores = NULL,
  chunkSize = 1000,
  verbose = TRUE
)
```

## Arguments

- complexes:

  A named list of protein complexes.

- enrichments:

  A named list of enrichment results (from \`runComplexEnrichment\`).

- mode:

  Edge weight mode: "compositional", "functional", or "combined".
  Defaults to "combined".

- similarityMethod:

  The metric for similarity. Defaults to \*\*"jaccard"\*\*.
  \*\*Warning:\*\* Using "overlap" is discouraged as it tends to
  collapse diverse hierarchies into single blobs.

- alpha:

  Numeric (0-1). Weight given to Compositional Similarity in "combined"
  mode. Defaults to \*\*0.75\*\* (favoring physical structure).

- nCores:

  Number of cores for parallel processing.

- chunkSize:

  Size of processing chunks.

- verbose:

  Logical.

## Value

A \`tibble\` of network edges with \`weight\` representing the
calculated similarity.

## Details

\*\*Systems Biology Rationale:\*\* To generate a diverse and physically
meaningful landscape, this function defaults to prioritizing
\*\*Compositional Similarity\*\* (shared proteins) over functional
similarity.

\- \*\*Composition (Protein Identity)\*\* is treated as the "Physical
Truth". - \*\*Function (Enrichment)\*\* is treated as an attribute.

The \`alpha\` parameter controls this balance. The new default (0.75)
gives 75 (like "Cell Cycle") from artificially pulling distinct physical
complexes into a single cluster, while still allowing functionally
related complexes to drift closer than unrelated ones.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
