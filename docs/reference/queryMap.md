# Query a ComplexMap Object for Specific Information

Allows for targeted querying of a \`ComplexMap\` object to find nodes
(complexes) that match specific criteria, such as containing a
particular protein or belonging to a biological theme.

## Usage

``` r
queryMap(complexMapObject, query, type)
```

## Arguments

- complexMapObject:

  A \`ComplexMap\` object.

- query:

  A character string containing the search term.

- type:

  A character string: one of \`"protein"\`, \`"complex"\`, or
  \`"theme"\`.

## Value

A \`tibble\` containing the rows from the node table that match the
query. Returns an empty tibble if no matches are found.

## Details

This function supports three distinct modes of querying: - \`type =
"protein"\`: Searches for complexes containing a specific protein. -
\`type = "complex"\`: Retrieves a specific complex by its ID. - \`type =
"theme"\`: Finds all complexes belonging to a given theme. This requires
theme information to be added to the \`ComplexMap\` object first (see
the "Advanced Analysis" vignette for an example).

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Assume 'cm_obj' is a valid ComplexMap object
# uba1_complexes <- queryMap(cm_obj, query = "UBA1", type = "protein")
```
