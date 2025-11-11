# Read a GMT File from a Local Path

Parses a local GMT (Gene Matrix Transposed) file into the standard
named-list format.

## Usage

``` r
getGmtFromFile(filepath, verbose = TRUE)
```

## Arguments

- filepath:

  The path to the local .gmt file.

- verbose:

  A logical value indicating whether to print progress messages.

## Value

A named list where names are gene set names and values are character
vectors of genes.
