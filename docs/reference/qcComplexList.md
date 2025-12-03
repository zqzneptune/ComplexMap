# Perform Quality Control on a List of Protein Complexes

A diagnostic tool to assess complex size distribution and distinguish
between true redundancy (synonyms) and biological hierarchy (subsets).

## Usage

``` r
qcComplexList(
  complexList,
  redundancyThreshold = 0.9,
  subsetThreshold = 0.9,
  verbose = TRUE
)
```

## Arguments

- complexList:

  A list of character vectors (protein complexes).

- redundancyThreshold:

  Threshold for Jaccard similarity warning. Defaults to 0.9.

- subsetThreshold:

  Threshold for Simpson coefficient (subset) reporting. Defaults to 0.9.

- verbose:

  Logical.

## Value

Invisibly returns the input \`complexList\`.

## Details

\*\*Systems Biology Rationale:\*\* To maintain a diverse functional
landscape, it is critical to distinguish between two types of
similarity:

1\. \*\*Jaccard Similarity (Redundancy):\*\* Complex A and B are nearly
identical. \* \*Action:\* Merge these (using \`refineComplexList\`). 2.
\*\*Overlap/Simpson Similarity (Hierarchy):\*\* Complex A is a subset of
Complex B. \* \*Action:\* \*\*Keep these.\*\* Do not merge. These often
represent distinct biological states (e.g., Core Complex vs.
Holo-Complex).

This report calculates both metrics to guide your refinement strategy.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
