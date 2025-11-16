# Perform Quality Control on a List of Protein Complexes

This function performs a quality control analysis on a list of protein
complexes, delivering a user-friendly report with key statistics and
actionable warnings.

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

  A logical value indicating whether to print the QC report to the
  console. Defaults to \`TRUE\`.

## Value

Invisibly returns the original \`complexList\` object, allowing it to be
used in a pipeline.

## Details

The QC process involves three main steps, presented in a clear
report: 1. \*\*Basic Statistics:\*\* Reports the total number of
complexes and unique proteins.

2\. \*\*Size Distribution:\*\* Summarizes the number of proteins per
complex, highlighting the smallest, median, and largest complexes. It
will issue an in-report warning if complexes have fewer than 3 members.

3\. \*\*Redundancy Analysis:\*\* Calculates the Jaccard similarity for
all pairs of complexes. It reports the median and maximum similarity and
issues an in-report warning if any pairs exceed the
\`redundancyThreshold\`.

This function is intended as a \*\*diagnostic tool\*\*. It uses the
Jaccard index for its clear and intuitive interpretation of redundancy.
For actively merging complexes, see the more advanced
\`refineComplexList()\` function.

## See also

\`refineComplexList()\` for a function to actively merge redundant
complexes.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Create a sample list of protein complexes
complex1 <- c("A", "B", "C", "D")
complex2 <- c("A", "B", "C", "E") # Highly redundant with complex1
complex3 <- c("F", "G", "H")
complex4 <- c("I", "J")          # Small complex
sampleList <- list(
  C1 = complex1, C2 = complex2, C3 = complex3, C4 = complex4
)

# Run the quality control analysis to see the new report format
qcComplexList(sampleList)
#> 
#> --- Quality Control Report for Complex List ---
#> 
#> [1] Basic Statistics
#>   ✓ Total Complexes: 4
#>   ✓ Unique Proteins: 10
#> 
#> [2] Complex Size Distribution
#>   - Smallest Complex: 2 proteins
#>   - Median Complex:   4 proteins
#>   - Largest Complex:  4 proteins
#>   ! WARNING: Found 1 complex(es) with fewer than 3 members.
#>     -> These are often considered too small for robust analysis.
#> 
#> [3] Redundancy Analysis (Jaccard Index)
#>   - Median Similarity: 0.000
#>   - Max Similarity:    0.600
#>   ✓ No highly redundant pairs detected (>= 0.80).
#> 
#> --- QC Complete ---
```
