# Benchmark Complex Refinement Parameters

Evaluates how different merging thresholds affect the quality AND
diversity of the resulting complex map.

## Usage

``` r
benchmarkParameters(
  complexList,
  referenceComplexes,
  threshold_range = seq(0.5, 1, by = 0.05),
  similarityMethod = "jaccard",
  ...
)
```

## Arguments

- complexList:

  A list of predicted protein complexes to be refined.

- referenceComplexes:

  A list of reference (gold standard) complexes.

- threshold_range:

  Numeric vector of thresholds to test. Defaults to \`seq(0.5, 1.0, by =
  0.05)\`.

- similarityMethod:

  The metric used for merging. Defaults to \*\*"jaccard"\*\*.

- ...:

  Additional arguments passed to \`refineComplexList\` (e.g., minSize).

## Value

A \`tibble\` with columns: - \`mergeThreshold\`: The threshold tested. -
\`NumComplexes\`: Number of complexes remaining (Diversity metric). -
\`PPV\`: Positive Predictive Value (Precision). - \`Sn\`: Sensitivity
(Recall). - \`F1\`: Harmonic mean of PPV and Sn. - \`Acc\`: Accuracy
(Geometric mean of PPV and Sn). - \`MMR\`: Maximum Matching Ratio.

## Details

\*\*Systems Biology Rationale:\*\* When optimizing the
\`mergeThreshold\`, there is a trade-off between removing redundancy
(increasing PPV) and losing specific variants (decreasing Sensitivity
and Diversity).

This function returns the standard metrics (PPV, Sn, Acc, MMR) but also
adds \*\*NumComplexes\*\*. - \*\*Collapse Warning:\*\* If
\`NumComplexes\` drops precipitously between two thresholds, you have
likely crossed the point where distinct biological variants are being
merged into generic blobs. - \*\*Recommendation:\*\* Choose the
threshold that maximizes F1/MMR while maintaining a \`NumComplexes\`
count close to your expected biological complexity.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
