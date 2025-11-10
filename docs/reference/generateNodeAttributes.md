# Generate Node Attributes for a Complex Network

Creates a detailed attribute table for each complex, suitable for
network visualization. Attributes include protein count, a primary
functional domain, top enriched functions, and a unique color
representing the complex's functional profile.

## Usage

``` r
generateNodeAttributes(
  complexes,
  enrichments,
  similarityMethod = "jaccard",
  verbose = TRUE
)
```

## Arguments

- complexes:

  A named list of protein complexes.

- enrichments:

  A named list of enrichment results, typically from
  \`runComplexEnrichment\`. Each element should be a data frame with at
  least \`ID\`, \`Description\`, and \`p.adjust\` columns.

- similarityMethod:

  The distance/similarity method passed to \`philentropy::distance\` for
  clustering terms. Defaults to "jaccard".

- verbose:

  A logical value indicating whether to print progress messages.
  Defaults to \`TRUE\`.

## Value

A \`tibble\` where each row corresponds to a complex. The columns
include: \`complexId\`, \`proteinCount\`, \`proteins\`,
\`primaryFunctionalDomain\`, \`topEnrichedFunctions\`, \`colorHex\`, and
\`sizeMapping\`.

## Details

This function performs several steps to generate rich node
attributes: 1. It aggregates all enriched terms from the input
\`enrichments\` list. 2. A term-complex matrix is built, and terms are
clustered based on their co-occurrence in complexes using the specified
similarity metric. This groups related functional terms into "functional
domains". 3. A unique color is assigned to each functional domain using
a qualitative palette from \`RColorBrewer\`. 4. For each complex, it
determines the primary functional domain based on the most significant
enriched term. 5. A unique "blended" color is calculated for each
complex by mixing the colors of its associated domains, weighted by the
significance (-log10 p-value) of the enriched terms. 6. Basic attributes
like protein count and a list of proteins are also included.

Complexes with no significant enrichments are assigned a default
"Unenriched" domain and a grey color.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data ---
complexes <- list(
  Cpx1 = c("A", "B", "C"),
  Cpx2 = c("C", "D", "E"),
  Cpx3 = c("F", "G") # Unenriched
)
enrichments <- list(
  Cpx1 = data.frame(ID = "GO:1", Description = "Term A", p.adjust = 0.01),
  Cpx2 = data.frame(ID = "GO:2", Description = "Term B", p.adjust = 0.02)
)

# --- Generate Node Attributes ---
nodeAttrs <- generateNodeAttributes(complexes, enrichments)
#> Generating core node attributes (function and color)...
#>     -> Clustering terms using 'jaccard' similarity.
#> Metric: 'jaccard' with unit: 'log'; comparing: 2 vectors
print(nodeAttrs)
#> # A tibble: 3 × 7
#>   complexId proteinCount proteins primaryFunctionalDomain topEnrichedFunctions
#>   <chr>            <int> <chr>    <chr>                   <chr>               
#> 1 Cpx1                 3 A,B,C    Term A                  Term A              
#> 2 Cpx2                 3 C,D,E    Term B                  Term B              
#> 3 Cpx3                 2 F,G      Unenriched              NA                  
#> # ℹ 2 more variables: colorHex <chr>, sizeMapping <dbl>
```
