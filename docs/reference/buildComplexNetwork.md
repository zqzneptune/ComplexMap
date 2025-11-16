# Build a Complex-Complex Interaction Network

Constructs a network of complexes where edges represent similarity. The
similarity can be based on shared proteins (compositional), shared
functional annotations (functional), or a weighted combination of both.

## Usage

``` r
buildComplexNetwork(
  complexes,
  enrichments,
  mode = "combined",
  similarityMethod = "jaccard",
  alpha = 0.5,
  nCores = NULL,
  chunkSize = 1000,
  verbose = TRUE
)
```

## Arguments

- complexes:

  A named list of protein complexes.

- enrichments:

  A named list of enrichment results, corresponding to the
  \`complexes\`. Typically the output of \`runComplexEnrichment\`.

- mode:

  A character string specifying how to calculate the final edge weight.
  One of "functional", "compositional", or "combined".

- similarityMethod:

  The metric used for both compositional and functional similarity. One
  of "jaccard", "overlap", or "dice".

- alpha:

  A numeric value (0-1) used only in "combined" mode to weigh the
  compositional similarity score.

- nCores:

  The number of CPU cores for parallel processing. Defaults to one less
  than available.

- chunkSize:

  The number of complex pairs to process in each parallel chunk.

- verbose:

  A logical value indicating whether to show progress messages and a
  progress bar.

## Value

A \`tibble\` representing the network edges. Each row includes the
source and target_complex_id complexes, their similarity scores, shared
component counts, the final calculated \`weight\`, and the
\`similarity_mode\`.

## Details

This function calculates all pairwise similarities between complexes in
a list. It is highly optimized for large datasets by using chunking and
parallel processing via the \`future\` framework.

The final edge weight is determined by the \`mode\`:

\- \`"compositional"\`: Uses only the protein similarity score.

\- \`"functional"\`: Uses only the functional annotation similarity
score.

\- \`"combined"\`: Uses a weighted average: \`alpha \* compositional +
(1 - alpha) \* functional\`.

If the \`progressr\` package is installed, a progress bar will be
displayed during the parallel computation when \`verbose = TRUE\`.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data (from previous examples) ---
complexes <- list(
  Cpx1 = c("A", "B", "C", "D"),
  Cpx2 = c("A", "B", "C", "E"), # similar to Cpx1
  Cpx3 = c("F", "G", "H")
)
enrichments <- list(
  Cpx1 = data.frame(ID = c("GO:1", "GO:2")),
  Cpx2 = data.frame(ID = c("GO:1", "GO:3")), # functionally similar to Cpx1
  Cpx3 = data.frame(ID = c("GO:4"))
)

# --- Build Network (using 2 cores for the example) ---
network <- buildComplexNetwork(
  complexes, enrichments, mode = "combined", nCores = 2
)
#> Building complex network using 'jaccard' similarity...
#> Using 2 cores for parallel processing.
#> Processing 3 complex pairs...
#> Split into 1 chunks of up to 1000 pairs each.
#> Combining results from chunks...
#> Calculating final weights and filtering...
#> Network construction complete: 1 edges retained.
print(network)
#> # A tibble: 1 × 8
#>   source_complex_id target_complex_id compSim funcSim sharedProt sharedFunc
#>   <chr>             <chr>               <dbl>   <dbl>      <int>      <int>
#> 1 Cpx1              Cpx2                  0.6   0.333          3          1
#> # ℹ 2 more variables: weight <dbl>, similarity_mode <chr>
```
