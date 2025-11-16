# 1. Preparing Your Input Data

## The ComplexMap Philosophy: Identifier Matching

The `ComplexMap` package was designed to be powerful, flexible, and
**species-agnostic**. It does not contain any hard-coded assumptions
about *Homo sapiens* or any other model organism.

The entire workflow depends on one simple principle: **the identifiers
used in your input lists must match**. This means the package’s
functions can support any organism and any identifier type (e.g., Gene
Symbols, Entrez IDs, Ensembl IDs, UniProt IDs), as long as they are
consistent.

There are three key inputs where identifiers must be consistent:

1.  **Your Complex List**: The list of protein/gene members for each
    complex you want to analyze.
2.  **Your Functional Gene Sets (GMT)**: The database used for
    enrichment analysis.
3.  **Your Reference Complex List**: A “gold standard” list used for
    benchmarking with
    [`ComplexMap::evaluateComplexes()`](https://zqzneptune.github.io/ComplexMap/reference/evaluateComplexes.md).
    (**This is only required for benchmarking, not for functional
    analysis.**)

For example, if your complex list uses **Gene Symbols**, your GMT file
must also use **Gene Symbols**. If your complex list uses **Entrez
IDs**, your GMT must use **Entrez IDs**.

This vignette demonstrates how to prepare the functional gene sets (GMT)
from various sources and, crucially, how to handle and convert
identifiers to ensure they match your input data.

``` r
# For this tutorial, we will assume our input complex list uses Gene Symbols.
myComplexes <- list(
  CPLX1 = c("POLR2A", "POLR2B", "POLR2C"),
  CPLX2 = c("CDK1", "CCNB1", "CCNB2")
)
```

## Preparing Functional Gene Sets (GMT)

`ComplexMap` provides several helper functions to obtain GMT files. To
align with the package’s philosophy of clarity and avoiding namespace
conflicts, we will call all functions explicitly using the
`package::function()` syntax (e.g.,
[`ComplexMap::getGmtFromFile()`](https://zqzneptune.github.io/ComplexMap/reference/getGmtFromFile.md))
instead of using [`library()`](https://rdrr.io/r/base/library.html).

Let’s explore each method, paying close attention to the identifier type
it returns.

### Method 1: From a User-Provided Local File

This is the most direct method. If you have your own GMT file, you can
load it with
[`ComplexMap::getGmtFromFile()`](https://zqzneptune.github.io/ComplexMap/reference/getGmtFromFile.md).

First, get the path to the example GMT file included with the package.
This file uses Gene Symbols.

``` r
gmtPath <- ComplexMap::getExampleGmt()
```

Load the GMT from the file path.

``` r
gmtFromFile <- ComplexMap::getGmtFromFile(gmtPath, verbose = FALSE)
```

Let’s inspect the identifiers:

``` r
# First 5 genes in the first gene set:
utils::head(gmtFromFile[[1]], 5)
#> [1] "ATF2"  "CHUK"  "IFNG"  "IKBKB" "IL2"
```

**Identifier Match:** The example GMT file uses **Gene Symbols**. Since
our hypothetical `myComplexes` list also uses Gene Symbols, these are
directly compatible and ready for analysis.

### Method 2: From the Molecular Signatures Database (MSigDB)

The `msigdbr` package provides a powerful and up-to-date interface to
the MSigDB collections. Our
[`ComplexMap::getMsigdbGmt()`](https://zqzneptune.github.io/ComplexMap/reference/getMsigdbGmt.md)
function simplifies this process.

``` r
# Fetch the Hallmark gene sets for Human
# This requires the `msigdbr` package
if (requireNamespace("msigdbr", quietly = TRUE)) {
  h_gmt <- ComplexMap::getMsigdbGmt(species = "Homo sapiens", collection = "H")

  # Inspect the identifiers
  utils::head(h_gmt[[1]], 5)
}
#> Fetching MSigDB sets (Species: Homo sapiens, Cat: H)
#> [1] "ABCA1" "ABCB8" "ACAA2" "ACADL" "ACADM"
```

**Identifier Match:** By default, `msigdbr` also returns **Gene
Symbols**. This is directly compatible with our `myComplexes` list.

### Method 3: From Gene Ontology (GO) via Bioconductor

Using official Bioconductor annotation packages is a highly reproducible
way to get functional annotations. These databases, however, typically
use stable database identifiers, not gene symbols.

``` r
# This requires an organism annotation package, e.g., org.Hs.eg.db for human
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  # Fetch Biological Process (BP) terms
  # We pass the database object directly using the `::` operator
  goGmt <- ComplexMap::getGoGmt(speciesDb = org.Hs.eg.db::org.Hs.eg.db, 
                                ontology = "BP",
                                verbose = FALSE)

  # Inspect the identifiers
  utils::head(goGmt[[1]], 5)
}
#> 
#> 'select()' returned 1:many mapping between keys and columns
#> 
#> 'select()' returned 1:1 mapping between keys and columns
#> [1] "1743" "2805" "2806" "3417" "3418"
```

**Identifier Mismatch!** The `getGoGmt` function returns a list where
the genes are **Entrez IDs** (e.g., “5594”, “5595”). These will **not
match** the Gene Symbols in our `myComplexes` list (e.g., “POLR2A”).

**Solution: Convert Your Complex List Identifiers**

The recommended approach is to convert your input complex identifiers to
match the stable IDs from the annotation database. The
[`AnnotationDbi::mapIds`](https://rdrr.io/pkg/AnnotationDbi/man/AnnotationDb-class.html)
function is perfect for this.

``` r
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  # Get all unique symbols from our complex list
  allSymbols <- unique(unlist(myComplexes))

  # Map symbols to Entrez IDs, taking the first ID if multiple exist
  symbolToEntrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = allSymbols,
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
  )
  
  # Remove any symbols that could not be mapped
  symbolToEntrez <- symbolToEntrez[!is.na(symbolToEntrez)]

  # Now, create a new complex list with Entrez IDs
  myComplexesEntrez <- lapply(myComplexes, function(complex) {
    # Look up the Entrez ID for each symbol
    entrez_ids <- symbolToEntrez[complex]
    # Return only the successfully mapped IDs, removing any NAs
    unname(entrez_ids[!is.na(entrez_ids)])
  })

  # Inspect the result
  print(myComplexesEntrez)
}
#> 'select()' returned 1:1 mapping between keys and columns
#> $CPLX1
#> [1] "5430" "5431" "5432"
#> 
#> $CPLX2
#> [1] "983"  "891"  "9133"
```

Now, the `myComplexesEntrez` list is directly compatible with the
`goGmt` generated from Bioconductor.

### Method 4: From Reactome Pathways via Bioconductor

Similarly, the `reactome.db` package provides pathway annotations, which
also use Entrez IDs.

``` r
if (requireNamespace("reactome.db", quietly = TRUE) && 
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      
  reactomeGmt <- ComplexMap::getReactomeGmt(
      speciesDb = org.Hs.eg.db::org.Hs.eg.db,
      verbose = FALSE
  )

  # Inspect the identifiers
  utils::head(reactomeGmt[[1]], 5)
}
#> [1] "1"     "10019" "10112" "10125" "10125"
```

**Identifier Mismatch!** Like the GO example, `reactome.db` provides
**Entrez IDs**.

**Solution:** The solution is the same as for Gene Ontology. You would
use the `myComplexesEntrez` list that we created in the previous step,
as its identifiers will match the identifiers in the `reactomeGmt`.

## Conclusion

This vignette has demonstrated the core philosophy of `ComplexMap`:
flexibility through identifier consistency. By understanding the
identifier types returned by different sources and knowing how to
convert your own data to match, you can apply the `ComplexMap` workflow
to virtually any organism.

Always check your identifiers before running
[`ComplexMap::createComplexMap()`](https://zqzneptune.github.io/ComplexMap/reference/createComplexMap.md)
to ensure a smooth and successful analysis. The explicit
`package::function()` syntax demonstrated here is a recommended best
practice for writing clean, reproducible, and error-free analysis
scripts.
