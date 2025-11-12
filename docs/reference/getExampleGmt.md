# Get the Path to an Example GMT File

Provides the full system path to an example gene set (GMT) file included
with the ComplexMap package.

## Usage

``` r
getExampleGmt()
```

## Value

A character string containing the full path to the example GMT file.

## Details

The included file is the BioCarta gene set collection from the Molecular
Signatures Database (MSigDB v2025.1). This function makes it easy to
access the file for use in examples and vignettes.

## Examples

``` r
# Get the path
gmtPath <- getExampleGmt()

# You can then read the file using the path
if (file.exists(gmtPath)) {
  exampleGmt <- getGmtFromFile(gmtPath)
}
#> Fetching gene sets from local file: /private/var/folders/x8/ngpdbcfd3pj5sf3r9kkxgs740000gn/T/RtmpojqfcE/temp_libpath40eca16c79e/ComplexMap/extdata/c2.cp.biocarta.v2025.1.Hs.symbols.gmt
```
