# ComplexMap <img src="man/figures/ComplexMap_logo.png" alt="ComplexMap hex logo" align="right" width ="139" />

**ComplexMap** provides an end-to-end workflow for **quality control**, **refinement**, **functional enrichment**, and **network-based analysis** of protein complex datasets.  

With this package, you can:  

- Quantitatively evaluate predicted complexes against a reference set  

- Perform functional enrichment analysis  

- Generate **static** and **interactive visualizations** to explore the functional landscape of protein complex maps  

Designed to take you from a **raw list of putative protein complexes** to a **functionally annotated, publication-ready network visualization**, ComplexMap streamlines complex analysis for researchers.

---

## **Installation**

You can install the development version of **ComplexMap** from GitHub.

### **1. Install Dependencies**

Before installing ComplexMap, make sure you have the required packages from **CRAN** and **Bioconductor**:

```r
# CRAN dependencies
install.packages(c(
  "dplyr", "tibble", "magrittr", "rlang", "stringr", "Matrix", 
  "igraph", "ggraph", "tidygraph", "ggplot2", "ggrepel", 
  "RColorBrewer", "scales", "clue", "future", "future.apply", 
  "philentropy", "visNetwork", "devtools"
))

# Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "msigdbr", "AnnotationDbi", "GO.db", "reactome.db", 
  "org.Hs.eg.db" # Example organism database
))
```

### **2. Install ComplexMap**

Once dependencies are installed, install ComplexMap from GitHub:

```r
devtools::install_github("zqzneptune/ComplexMap")
```

---

## **Quick Start Example**

Here’s a minimal example of the core workflow: from a list of complexes to a functional network map.

```r
library(ComplexMap)

# 1. Load example data
data("demoComplexes")
gmtPath <- getExampleGmt()
gmt <- getGmtFromFile(gmtPath)

# 2. Refine complexes
refinedComplexes <- refineComplexList(demoComplexes, mergeThreshold = 0.75, verbose = FALSE)

# 3. Functional enrichment
enrichments <- runComplexEnrichment(refinedComplexes, gmt, verbose = FALSE)

# 4. Build network
networkEdges <- buildComplexNetwork(refinedComplexes, enrichments, verbose = FALSE)

# 5. Generate node attributes and layout
nodeAttributes <- generateNodeAttributes(refinedComplexes, enrichments, verbose = FALSE)
mapLayout <- computeMapTopology(nodeAttributes, networkEdges, verbose = FALSE)

# 6. Visualize the final map
visualizeMapWithLegend(mapLayout, networkEdges, verbose = FALSE)
```

---

## **Learn More**

For a **full walkthrough**, detailed function documentation, and tutorials, visit:  
👉 https://zqzneptune.github.io/ComplexMap/
