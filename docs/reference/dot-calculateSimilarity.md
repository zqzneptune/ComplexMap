# Calculate a Similarity Score Between Two Sets

An internal helper function that computes one of three common similarity
metrics: Jaccard, Overlap, or Dice.

## Usage

``` r
.calculateSimilarity(set1, set2, method = "jaccard")
```

## Arguments

- set1:

  A character vector representing the first set.

- set2:

  A character vector representing the second set.

- method:

  The similarity metric to use. One of "jaccard", "overlap", or "dice".

## Value

A numeric similarity score between 0 and 1. Returns 0 if either set is
empty or if there is no intersection.
