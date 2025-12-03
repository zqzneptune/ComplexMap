# Export a ComplexMap Network for External Tools (Cytoscape)

Exports the \`ComplexMap\` nodes and edges to compatible file formats,
preserving systems biology attributes like Specificity Scores and
Functional Domains.

## Usage

``` r
exportNetwork(
  complexMapObject,
  filePrefix,
  format = "cytoscape",
  verbose = TRUE
)
```

## Arguments

- complexMapObject:

  A \`ComplexMap\` object.

- filePrefix:

  Character string for output filenames (e.g., "results/my_map").

- format:

  Output format. Currently supports "cytoscape" (TSV).

- verbose:

  Logical.

## Value

None (Writes files to disk).

## Details

\*\*Systems Biology Rationale:\*\* To verify the landscape in external
tools (e.g., Cytoscape), this function ensures that the specific
calculated attributes are correctly formatted:

\- \*\*\`score\`\*\*: The Specificity Score (Size/Color mapping). -
\*\*\`primaryFunctionalDomain\`\*\*: The Specific Label. -
\*\*\`colorHex\`\*\*: The calculated blended color.

It performs a "Sanitization" step to convert any R-specific list columns
into semi-colon separated strings and fills \`NA\` values to ensure safe
import.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>
