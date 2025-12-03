# Print Method for a ComplexMap Object

Provides a "Systems Biology Dashboard" summary of the ComplexMap object.

## Usage

``` r
# S3 method for class 'ComplexMap'
print(x, ...)
```

## Arguments

- x:

  A \`ComplexMap\` object.

- ...:

  Additional arguments.

## Value

Invisibly returns the object.

## Details

This summary is designed to give immediate feedback on the
\*\*Functional Diversity\*\* of the landscape:

\- \*\*Functional Diversity:\*\* The count of distinct functional labels
(colors). A high number (e.g., \>15) indicates a rich, specific
landscape. A low number indicates over-clustering or generic
annotation. - \*\*Annotation Coverage:\*\* The percentage of complexes
that could be assigned a function. High coverage with high diversity is
the ideal state.
