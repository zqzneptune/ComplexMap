# Read a GMT File from a Local Path

Parses a local GMT (Gene Matrix Transposed) file into the standard
named-list format.

## Usage

``` r
getGmtFromFile(filepath)
```

## Arguments

- filepath:

  The path to the local .gmt file.

## Value

A named list where names are gene set names and values are character
vectors of genes.
