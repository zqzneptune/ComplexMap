# Refine a List of Protein Complexes by Size and Redundancy

This function refines a list of protein complexes by first filtering
them based on size and then merging highly redundant complexes based on
a Jaccard similarity threshold.

## Usage

``` r
refineComplexList(
  complexList,
  minSize = 3,
  maxSize = 500,
  mergeThreshold = 0.9,
  verbose = TRUE
)
```

## Arguments

- complexList:

  A named list where each element is a character vector of protein
  identifiers representing a complex.

- minSize:

  An integer specifying the minimum number of proteins a complex must
  have to be retained. Defaults to 3.

- maxSize:

  An integer specifying the maximum number of proteins a complex can
  have to be retained. Defaults to 500.

- mergeThreshold:

  A numeric value (0-1) for the Jaccard similarity. Complexes with a
  score \>= this value will be merged. Defaults to 0.9.

- verbose:

  A logical value indicating whether to print progress messages.
  Defaults to \`TRUE\`.

## Value

A refined and renamed named list of protein complexes.

## Details

The refinement process consists of two main stages: 1. \*\*Size
Filtering:\*\* Complexes that are smaller than \`minSize\` or larger
than \`maxSize\` are removed from the list. 2. \*\*Redundancy
Merging:\*\* A Jaccard similarity matrix is calculated for all remaining
pairs of complexes. A union-find algorithm is then used to identify
clusters (or "merge groups") of complexes where all members are
connected by a similarity score greater than or equal to the
\`mergeThreshold\`. These groups are then merged into single complexes.

Finally, all complexes in the refined list are renamed to a standardized
format ("CpxMap_0001", "CpxMap_0002", etc.).

## See also

\`qcComplexList()\` for quality control analysis of a complex list.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Create a sample list of protein complexes
complex1 <- c("A", "B", "C", "D")
complex2 <- c("A", "B", "C", "E") # Highly redundant with complex1
complex3 <- c("F", "G", "H")
complex4 <- c("I", "J")          # Too small, will be filtered
sampleList <- list(
  C1 = complex1, C2 = complex2, C3 = complex3, C4 = complex4
)

# Refine the list using a high similarity threshold
refinedList <- refineComplexList(sampleList, mergeThreshold = 0.75)
#> 
#> --- Refining Input Complex List ---
#> Filtered 1 complexes by size. Retaining 3.
#> Identifying merge groups with Jaccard >= 0.75...
#> Found 0 merge groups. Merging 0 complexes into 3.
#> Merging complete. Final list has 3 complexes.
#> 
#> --- Refinement Complete ---
```
