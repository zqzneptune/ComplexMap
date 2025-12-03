# Query a ComplexMap Object

Targeted querying of complexes, proteins, or themes.

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

  "protein", "complex", or "theme".

## Value

A tibble of matching nodes.

## Details

\- \`type="protein"\`: Regex search for protein members. -
\`type="complex"\`: Exact match for Complex ID. - \`type="theme"\`:
Exact match for Theme Label. Requires \`summarizeThemes()\` to be run on
the object first with \`add_to_object = TRUE\`.
