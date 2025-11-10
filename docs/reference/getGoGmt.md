# Fetch Gene Ontology (GO) Gene Sets

Fetches GO terms and their associated genes from a Bioconductor
AnnotationDb package.

## Usage

``` r
getGoGmt(speciesDb, ontology = "BP", minGmtSize = 10, maxGmtSize = 500)
```

## Arguments

- speciesDb:

  An AnnotationDb object for the target species (e.g.,
  \`org.Hs.eg.db\`).

- ontology:

  The GO ontology to fetch. One of "BP", "MF", or "CC".

- minGmtSize:

  The minimum number of genes for a GO term to be included.

- maxGmtSize:

  The maximum number of genes for a GO term to be included.

## Value

A named list where names are GO terms and values are character vectors
of Entrez IDs.
