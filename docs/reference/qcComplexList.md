# Perform Quality Control on a List of Protein Complexes

This function performs a quality control analysis on a list of protein
complexes. It provides basic statistics, analyzes the distribution of
complex sizes, and calculates pairwise redundancy using the Jaccard
similarity index.

## Usage

``` r
qcComplexList(complexList, redundancyThreshold = 0.8, verbose = TRUE)
```

## Arguments

- complexList:

  A list where each element is a character vector of protein identifiers
  representing a complex.

- redundancyThreshold:

  A numeric value between 0 and 1. A warning will be issued for any pair
  of complexes with a Jaccard similarity score greater than or equal to
  this threshold. Defaults to 0.8.

- verbose:

  A logical value indicating whether to print progress messages and
  summaries to the console. Defaults to \`TRUE\`.

## Value

Invisibly returns the original \`complexList\` object, allowing it to be
used in a pipeline.

## Details

The QC process involves three main steps: 1. \*\*Basic Statistics:\*\*
Reports the total number of complexes and unique proteins in the list.
2. \*\*Size Distribution:\*\* Provides a summary of the number of
proteins per complex and warns if complexes have fewer than three
members. 3. \*\*Redundancy Analysis:\*\* A sparse binary membership
matrix is constructed to efficiently calculate the Jaccard similarity
for all unique pairs of complexes. A warning is issued if a significant
portion of pairs exceeds the specified \`redundancyThreshold\`.

## See also

\`refine_complex_list()\` for a method to merge redundant complexes.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Create a sample list of protein complexes
complex1 <- c("A", "B", "C", "D")
complex2 <- c("A", "B", "C", "E") # Highly redundant with complex1
complex3 <- c("F", "G", "H")
complex4 <- c("I", "J")          # Small complex
complex5 <- c("X", "Y", "Z")
sampleList <- list(
  C1 = complex1, C2 = complex2, C3 = complex3, C4 = complex4, C5 = complex5
)

# Run the quality control analysis
qcComplexList(sampleList)
#> 
#> --- Running Quality Control on Complex List ---
#> 
#> [1] Basic Statistics:
#>     - Total number of complexes: 5
#>     - Total number of unique proteins: 13
#> 
#> [2] Complex Size Distribution:
#>        Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>         2.0     3.0     3.0     3.2     4.0     4.0 
#> Warning: 1 complexes have fewer than 3 members.
#> 
#> [3] Redundancy Analysis:
#>     - Distribution of Jaccard similarity scores:
#>        Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>        0.00    0.00    0.00    0.06    0.00    0.60 
#>     - No highly redundant complex pairs detected.
#> 
#> --- QC Complete ---
```
