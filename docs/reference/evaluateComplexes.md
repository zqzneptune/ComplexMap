# Evaluate Predicted Protein Complexes Against a Reference Set

Calculates four standard metrics for evaluating protein complex
predictions: Positive Predictive Value (PPV), Sensitivity (Sn), Accuracy
(Acc), and the Maximum Matching Ratio (MMR).

## Usage

``` r
evaluateComplexes(
  predictedComplexes,
  referenceComplexes,
  nCores = NULL,
  verbose = TRUE
)
```

## Arguments

- predictedComplexes:

  A list of predicted protein complexes.

- referenceComplexes:

  A list of reference (gold standard) complexes.

- nCores:

  The number of CPU cores to use for parallel computation. Defaults to
  one less than the total number of detected cores.

- verbose:

  A logical value indicating whether to print progress messages.
  Defaults to \`TRUE\`.

## Value

A named list containing four numeric values: \`PPV\`, \`Sn\`, \`Acc\`,
and \`MMR\`. Returns \`NA\` for all metrics if either input list is
empty.

## Details

This function is optimized for speed by calculating a shared
intersection matrix between predicted and reference complexes in
parallel. This matrix is then used as the basis for all four metric
calculations.

\- \*\*PPV, Sn, and Acc\*\* are calculated based on the confusion matrix
between predicted and reference complexes, as described in the
literature (e.g., Zhang et al., 2012).

\- \*\*MMR\*\* is calculated by first deriving an overlap score matrix,
where the score for a predicted complex (P) and a reference complex (R)
is \`\|P ∩ R\|² / (\|P\| \* \|R\|)\`. The \[Hungarian
algorithm\](https://en.wikipedia.org/wiki/Hungarian_algorithm) is then
used to solve the maximum weight bipartite matching problem. This
approach is based on the method described by Nepusz et al. (2012).

The parallel computation uses \`parallel::mclapply\`, which is not
available on Windows. On Windows, the calculation will run sequentially.

## References

Nepusz, T., Yu, H. & Paccanaro, A. (2012). Detecting overlapping protein
complexes in protein-protein interaction networks. \*Nature Methods\*,
9, 471–472. [doi:10.1038/nmeth.1938](https://doi.org/10.1038/nmeth.1938)

Zhang XF, Dai DQ, Ou-Yang L, Wu MY (2012). Exploring Overlapping
Functional Units with Various Structure in Protein Interaction Networks.
\*PLOS ONE\*, 7(8): e43092.
[doi:10.1371/journal.pone.0043092](https://doi.org/10.1371/journal.pone.0043092)

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data ---
# Predicted complexes
pred1 <- c("A", "B", "C")
pred2 <- c("D", "E", "F")
pred3 <- c("A", "G", "H")
predicted <- list(P1 = pred1, P2 = pred2, P3 = pred3)

# Reference complexes (gold standard)
ref1 <- c("A", "B", "C", "X") # Good match for pred1
ref2 <- c("D", "E", "F")     # Perfect match for pred2
ref3 <- c("I", "J", "K")     # Unmatched complex
reference <- list(R1 = ref1, R2 = ref2, R3 = ref3)

# --- Evaluation ---
# Use 2 cores for the example
metrics <- evaluateComplexes(predicted, reference, nCores = 2)
#> 
#> --- Evaluating Complex Predictions ---
#> [1] Calculating intersections using 2 core(s)...
#> [2] Calculating PPV, Sensitivity, and Accuracy...
#> [3] Calculating Maximum Matching Ratio (MMR)...
#> 
#> --- Evaluation Complete ---
print(metrics)
#> $PPV
#> [1] 1
#> 
#> $Sn
#> [1] 0.6
#> 
#> $Acc
#> [1] 0.7745967
#> 
#> $MMR
#> [1] 0.5833333
#> 
```
