# Refine a List of Protein Complexes by Size and Redundancy

This function filters complexes by size and merges highly redundant
ones.

## Usage

``` r
refineComplexList(
  complexList,
  minSize = 3,
  maxSize = 500,
  mergeThreshold = 0.9,
  similarityMethod = "jaccard",
  verbose = TRUE
)
```

## Arguments

- complexList:

  A named list where each element is a character vector of protein
  identifiers representing a complex.

- minSize:

  Minimum complex size. Defaults to 3.

- maxSize:

  Maximum complex size. Defaults to 500.

- mergeThreshold:

  Numeric (0-1). Threshold for merging. Defaults to \*\*0.90\*\*
  (strict) to preserve distinct variants.

- similarityMethod:

  Metric for merging. Defaults to \*\*"jaccard"\*\*. Options: \*
  \`"jaccard"\`: (Default) \`Intersection / Union\`. Penalizes size
  differences. Best for preserving specific variants. \* \`"overlap"\`:
  \`Intersection / Min(A,B)\`. Merges subsets into parents. Use only if
  you want to aggressively reduce data. \* \`"matching_score"\`:
  \`Intersection^2 / (SizeA \* SizeB)\`. Geometric approach. \*
  \`"dice"\`: \`2 \* Intersection / (SizeA + SizeB)\`.

- verbose:

  Logical.

## Value

A list containing: \* \`refinedComplexes\`: The final list of merged
complexes. \* \`mergeMap\`: A tibble mapping original IDs to final IDs.

## Details

\*\*Systems Biology Rationale:\*\* To preserve the functional diversity
of the input landscape, this function defaults to a conservative
\*\*Jaccard\*\* similarity with a high threshold (0.90).

\* \*\*Trust the Input:\*\* We assume the upstream clustering method
identified biologically relevant variants (e.g., a complex with vs.
without a regulatory subunit). \* \*\*Minimal Merging:\*\* We only merge
complexes if they are effectively synonyms (nearly identical
composition).

We discourage using "overlap" or "matching_score" for this step, as they
tend to merge sub-complexes into parents, destroying the specific
functional signals you wish to visualize.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
