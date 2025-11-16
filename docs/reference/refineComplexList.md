# Refine a List of Protein Complexes by Size and Redundancy

This function refines a list of protein complexes by first filtering
them based on size and then merging highly redundant complexes based on
a similarity threshold.

## Usage

``` r
refineComplexList(
  complexList,
  minSize = 3,
  maxSize = 500,
  mergeThreshold = 0.9,
  similarityMethod = "matching_score",
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

  A numeric value (0-1). Complexes with a similarity score \>= this
  value will be merged. Defaults to 0.9.

- similarityMethod:

  A character string specifying the similarity metric for merging.
  Defaults to \`"matching_score"\`. The available options are:

  \`"matching_score"\`

  :   \`Intersection² / (\|A\| \* \|B\|)\`. Rewards large, shared cores.
      Aligns with the MMR evaluation metric.

  \`"simpson"\`

  :   \`Intersection / min(\|A\|, \|B\|)\`. Also known as the Overlap
      coefficient. excels at merging sub-complexes into larger parents.

  \`"jaccard"\`

  :   \`Intersection / Union\`. A classic, balanced metric that
      penalizes size differences.

  \`"dice"\`

  :   \`2 \* Intersection / (\|A\| + \|B\|)\`. Similar to Jaccard but
      generally less stringent.

- verbose:

  A logical value indicating whether to print progress messages.
  Defaults to \`TRUE\`.

## Value

A list containing two named elements:

- \`refinedComplexes\`:

  The final, renamed list of merged complexes.

- \`mergeMap\`:

  A \`tibble\` with columns \`originalId\` and \`finalId\`, mapping each
  original complex to its final standardized ID.

## Details

The refinement process consists of two main stages:

1\. \*\*Size Filtering:\*\* Complexes smaller than \`minSize\` or larger
than \`maxSize\` are removed.

2\. \*\*Redundancy Merging:\*\* A similarity matrix is calculated for
all remaining complex pairs using the chosen \`similarityMethod\`. A
union-find algorithm identifies clusters of complexes connected by a
similarity score \>= \`mergeThreshold\`, which are then merged.

Finally, all complexes in the refined list are renamed to a standardized
format ("CpxMap_0001", etc.), and a traceability map is generated.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Create a sample list of protein complexes
c1 <- c("A", "B", "C", "D", "E")
c2 <- c("A", "B", "C", "D", "F") # High similarity with c1
c3 <- c("A", "B", "C")          # Subset of c1 and c2
c4 <- c("X", "Y", "Z")
sampleList <- list(C1=c1, C2=c2, C3=c3, C4=c4)

# Refine using the default "matching_score"
refineComplexList(sampleList, mergeThreshold = 0.6)
#> 
#> --- Refining Input Complex List ---
#> Filtered 0 complexes by size. Retaining 4.
#> Identifying merge groups with matching_score >= 0.60...
#> Found 1 merge groups. Merging 2 complexes into 2.
#> Merging complete. Final list has 2 complexes.
#> 
#> --- Refinement Complete ---
#> $refinedComplexes
#> $refinedComplexes$CpxMap_0001
#> [1] "A" "B" "C" "D" "E" "F"
#> 
#> $refinedComplexes$CpxMap_0002
#> [1] "X" "Y" "Z"
#> 
#> 
#> $mergeMap
#> # A tibble: 4 × 2
#>   originalId finalId    
#>   <chr>      <chr>      
#> 1 C1         CpxMap_0001
#> 2 C2         CpxMap_0001
#> 3 C3         CpxMap_0001
#> 4 C4         CpxMap_0002
#> 

# Refine using the "simpson" method to merge the subset
refineComplexList(sampleList, mergeThreshold = 0.9, similarityMethod = "simpson")
#> 
#> --- Refining Input Complex List ---
#> Filtered 0 complexes by size. Retaining 4.
#> Identifying merge groups with simpson >= 0.90...
#> Found 1 merge groups. Merging 2 complexes into 2.
#> Merging complete. Final list has 2 complexes.
#> 
#> --- Refinement Complete ---
#> $refinedComplexes
#> $refinedComplexes$CpxMap_0001
#> [1] "A" "B" "C" "D" "E" "F"
#> 
#> $refinedComplexes$CpxMap_0002
#> [1] "X" "Y" "Z"
#> 
#> 
#> $mergeMap
#> # A tibble: 4 × 2
#>   originalId finalId    
#>   <chr>      <chr>      
#> 1 C1         CpxMap_0001
#> 2 C2         CpxMap_0001
#> 3 C3         CpxMap_0001
#> 4 C4         CpxMap_0002
#> 
```
