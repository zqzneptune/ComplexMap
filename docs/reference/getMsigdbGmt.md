# Fetch Gene Sets from MSigDB

A wrapper for \`msigdbr::msigdbr\` to fetch gene sets and format them as
a named list (GMT format).

## Usage

``` r
getMsigdbGmt(species = "Homo sapiens", collection = "H", verbose = TRUE)
```

## Arguments

- species:

  The scientific name for the species (e.g., "Homo sapiens").

- collection:

  The MSigDB collection code (e.g., "H" for hallmark).

- verbose:

  A logical value indicating whether to print progress messages.

## Value

A named list where names are gene set names and values are character
vectors of gene symbols.
