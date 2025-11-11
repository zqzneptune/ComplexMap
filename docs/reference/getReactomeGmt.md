# Fetch Reactome Pathway Gene Sets

Fetches Reactome pathways and their associated genes from the
\`reactome.db\` Bioconductor annotation package.

## Usage

``` r
getReactomeGmt(speciesDb, verbose = TRUE)
```

## Arguments

- speciesDb:

  An AnnotationDb object for the target species (e.g.,
  \`org.Hs.eg.db\`). This is used to filter pathways for the correct
  species.

- verbose:

  A logical value indicating whether to print progress messages.

## Value

A named list where names are Reactome pathway names and values are
character vectors of Entrez IDs.
