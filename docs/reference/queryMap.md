# Query a ComplexMap Object

Targeted querying of complexes or proteins.

## Usage

``` r
queryMap(complexMapObject, query, type)
```

## Arguments

- complexMapObject:

  A \`ComplexMap\` object.

- query:

  Search string.

- type:

  "protein" or "complex".

## Value

A tibble of matching nodes.

## Details

\- \`type="protein"\`: Regex search for protein members. -
\`type="complex"\`: Exact match for Complex ID.
