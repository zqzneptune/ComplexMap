# Export a ComplexMap Network for External Tools

This function exports the node and edge tables from a \`ComplexMap\`
object into a format suitable for widely used network visualization
software, such as Cytoscape.

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

  A \`ComplexMap\` object returned by \`createComplexMap()\`.

- filePrefix:

  A character string that will be used as the prefix for the output
  filenames. For example, a \`filePrefix\` of "my_map" will result in
  "my_map_nodes.tsv" and "my_map_edges.tsv".

- format:

  A character string specifying the output format. Currently, only
  \`"cytoscape"\` is supported.

- verbose:

  A logical value indicating whether to print a confirmation message
  upon successful export.

## Value

The function is called for its side effect of writing files to disk. It
does not return a value.

## Details

The function currently supports one format: - \`"cytoscape"\`: This
option generates two separate tab-separated value (.tsv) files. One file
contains the node attributes (\`\<filePrefix\>\_nodes.tsv\`) and the
other contains the edge list with its attributes
(\`\<filePrefix\>\_edges.tsv\`). These files can be directly imported
into Cytoscape's network and attribute tables.

The function uses \`utils::write.table\` for robust and
standard-compliant file writing.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# Assume 'cm_obj' is a valid ComplexMap object
# dir <- tempdir() # Use a temporary directory for the example
# exportNetwork(cm_obj, filePrefix = file.path(dir, "myNetwork"))
```
