# Preparing Input Data

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
    [`evaluateComplexes()`](https://zqzneptune.github.io/ComplexMap/reference/evaluateComplexes.md).
    (**This is only required for benchmarking, not for functional
    analysis.**)

For example, if your complex list uses **Gene Symbols**, your GMT file
must also use **Gene Symbols**. If your complex list uses **Entrez
IDs**, your GMT must use **Entrez IDs**.

This vignette demonstrates how to prepare the functional gene sets (GMT)
from various sources and, crucially, how to handle and convert
identifiers to ensure they match your input data.

``` r

library(ComplexMap)
library(dplyr)

# For this tutorial, we will assume our input complex list uses Gene Symbols.
myComplexes <- list(
  CPLX1 = c("POLR2A", "POLR2B", "POLR2C"),
  CPLX2 = c("CDK1", "CCNB1", "CCNB2")
)
```

## Preparing Functional Gene Sets (GMT)

`ComplexMap` provides several helper functions to obtain GMT files.
Let’s explore each one, paying close attention to the identifier type it
returns.

### Method 1: From a User-Provided Local File

This is the most direct method. If you have your own GMT file, you can
load it with
[`getGmtFromFile()`](https://zqzneptune.github.io/ComplexMap/reference/getGmtFromFile.md).

First, get the path to the example GMT file included with the package.
This file uses Gene Symbols.

``` r

gmtPath <- getExampleGmt()
```

Load the GMT from the file path

``` r

gmtFromFile <- getGmtFromFile(gmtPath)
```

    ## Fetching gene sets from local file: /private/var/folders/gb/q0_2jm654r9_t3r2hb11v3tm0000gn/T/Rtmpx4UtqW/temp_libpath4397f75503b/ComplexMap/extdata/c2.cp.biocarta.v2025.1.Hs.symbols.gmt

Let’s inspect the identifiers

``` r

# Name of the first gene set:
names(gmtFromFile)[1:5]
```

    ## [1] "BIOCARTA_41BB_PATHWAY"          "BIOCARTA_ACE2_PATHWAY"         
    ## [3] "BIOCARTA_ACETAMINOPHEN_PATHWAY" "BIOCARTA_ACH_PATHWAY"          
    ## [5] "BIOCARTA_ACTINY_PATHWAY"

``` r

# First 5 genes in that set:
gmtFromFile[1]
```

    ## $BIOCARTA_41BB_PATHWAY
    ##  [1] "ATF2"    "CHUK"    "IFNG"    "IKBKB"   "IL2"     "IL4"     "JUN"    
    ##  [8] "MAP3K1"  "MAP3K5"  "MAP4K5"  "MAPK14"  "MAPK8"   "NFKB1"   "NFKBIA" 
    ## [15] "RELA"    "TNFRSF9" "TNFSF9"  "TRAF2"

**Identifier Match:** The example GMT file uses **Gene Symbols**. Since
our hypothetical `myComplexes` list also uses Gene Symbols, these are
directly compatible and ready for analysis.

### Method 2: From the Molecular Signatures Database (MSigDB)

The `msigdbr` package provides a powerful and up-to-date interface to
the MSigDB collections. Our
[`getMsigdbGmt()`](https://zqzneptune.github.io/ComplexMap/reference/getMsigdbGmt.md)
function simplifies this process.

``` r

# Fetch the Hallmark gene sets for Human
# This requires the `msigdbr` package
if (requireNamespace("msigdbr", quietly = TRUE)) {
  h_gmt <- getMsigdbGmt(species = "Homo sapiens", collection = "H")

  # Inspect the identifiers
  h_gmt[1]
}
```

    ## Fetching MSigDB sets (Species: Homo sapiens, Cat: H)

    ## $HALLMARK_ADIPOGENESIS
    ##   [1] "ABCA1"    "ABCB8"    "ACAA2"    "ACADL"    "ACADM"    "ACADS"   
    ##   [7] "ACLY"     "ACO2"     "ACOX1"    "ADCY6"    "ADIG"     "ADIPOQ"  
    ##  [13] "ADIPOR2"  "AGPAT3"   "AIFM1"    "AK2"      "ALDH2"    "ALDOA"   
    ##  [19] "ANGPT1"   "ANGPTL4"  "APLP2"    "APOE"     "ARAF"     "ARL4A"   
    ##  [25] "ATL2"     "ATP1B3"   "ATP5PO"   "BAZ2A"    "BCKDHA"   "BCL2L13" 
    ##  [31] "BCL6"     "C3"       "CAT"      "CAVIN1"   "CAVIN2"   "CCNG2"   
    ##  [37] "CD151"    "CD302"    "CD36"     "CDKN2C"   "CHCHD10"  "CHUK"    
    ##  [43] "CIDEA"    "CMBL"     "CMPK1"    "COL15A1"  "COL4A1"   "COQ3"    
    ##  [49] "COQ5"     "COQ9"     "COX6A1"   "COX7B"    "COX8A"    "CPT2"    
    ##  [55] "CRAT"     "CS"       "CYC1"     "CYP4B1"   "DBT"      "DDT"     
    ##  [61] "DECR1"    "DGAT1"    "DHCR7"    "DHRS7"    "DHRS7B"   "DLAT"    
    ##  [67] "DLD"      "DNAJB9"   "DNAJC15"  "DRAM2"    "ECH1"     "ECHS1"   
    ##  [73] "ELMOD3"   "ELOVL6"   "ENPP2"    "EPHX2"    "ESRRA"    "ESYT1"   
    ##  [79] "ETFB"     "FABP4"    "FAH"      "FZD4"     "G3BP2"    "GADD45A" 
    ##  [85] "GBE1"     "GHITM"    "GPAM"     "GPAT4"    "GPD2"     "GPHN"    
    ##  [91] "GPX3"     "GPX4"     "GRPEL1"   "HADH"     "HIBCH"    "HSPB8"   
    ##  [97] "IDH1"     "IDH3A"    "IDH3G"    "IFNGR1"   "IMMT"     "ITGA7"   
    ## [103] "ITIH5"    "ITSN1"    "JAGN1"    "LAMA4"    "LEP"      "LIFR"    
    ## [109] "LIPE"     "LPCAT3"   "LPL"      "LTC4S"    "MAP4K3"   "MCCC1"   
    ## [115] "MDH2"     "ME1"      "MGLL"     "MGST3"    "MIGA2"    "MRAP"    
    ## [121] "MRPL15"   "MTARC2"   "MTCH2"    "MYLK"     "NABP1"    "NDUFA5"  
    ## [127] "NDUFAB1"  "NDUFB7"   "NDUFS3"   "NKIRAS1"  "NMT1"     "OMD"     
    ## [133] "ORM1"     "PDCD4"    "PEMT"     "PEX14"    "PFKFB3"   "PFKL"    
    ## [139] "PGM1"     "PHLDB1"   "PHYH"     "PIM3"     "PLIN2"    "POR"     
    ## [145] "PPARG"    "PPM1B"    "PPP1R15B" "PRDX3"    "PREB"     "PTCD3"   
    ## [151] "PTGER3"   "QDPR"     "RAB34"    "REEP5"    "REEP6"    "RETN"    
    ## [157] "RETSAT"   "RIOK3"    "RMDN3"    "RNF11"    "RREB1"    "RTN3"    
    ## [163] "SAMM50"   "SCARB1"   "SCP2"     "SDHB"     "SDHC"     "SLC19A1" 
    ## [169] "SLC1A5"   "SLC25A1"  "SLC25A10" "SLC27A1"  "SLC5A6"   "SLC66A3" 
    ## [175] "SNCG"     "SOD1"     "SORBS1"   "SOWAHC"   "SPARCL1"  "SQOR"    
    ## [181] "SSPN"     "STAT5A"   "STOM"     "SUCLG1"   "SULT1A1"  "TALDO1"  
    ## [187] "TANK"     "TKT"      "TOB1"     "TST"      "UBC"      "UBQLN1"  
    ## [193] "UCK1"     "UCP2"     "UQCR10"   "UQCR11"   "UQCRC1"   "UQCRQ"   
    ## [199] "VEGFB"    "YWHAG"

**Identifier Match:** By default, `msigdbr` also returns **Gene
Symbols**. This is directly compatible with our `myComplexes` list.

### Method 3: From Gene Ontology (GO) via Bioconductor

Using official Bioconductor annotation packages is a highly reproducible
way to get functional annotations. These databases, however, typically
use stable database identifiers, not gene symbols.

``` r

# This requires an organism annotation package, e.g., org.Hs.eg.db for human
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  # Fetch Biological Process (BP) terms
  goGmt <- getGoGmt(speciesDb = org.Hs.eg.db, ontology = "BP")

  # Inspect the identifiers
  goGmt[1]
}
```

    ## 

    ## Fetching GO (BP) gene sets from Bioconductor...

    ## 'select()' returned 1:many mapping between keys and columns

    ## 

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Retained 2942 GO terms (size 10-500) out of 12173

    ## $`2-oxoglutarate metabolic process`
    ##  [1] "1291"   "1738"   "1743"   "2805"   "2806"   "3417"   "3418"   "4967"  
    ##  [9] "5264"   "6898"   "51166"  "55753"  "56267"  "79944"  "84706"  "92259" 
    ## [17] "137872" "728294"

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
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  # Get all unique symbols from our complex list
  allSymbols <- unique(unlist(myComplexes))

  # Map symbols to Entrez IDs
  symbolToEntrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                           keys = allSymbols,
                                           keytype = "SYMBOL",
                                           column = "ENTREZID",
                                           multiVals = "first")
  
  # Remove any symbols that could not be mapped
  symbolToEntrez <- symbolToEntrez[!is.na(symbolToEntrez)]

  # Now, create a new complex list with Entrez IDs
  myComplexesEntrez <- lapply(myComplexes, function(complex) {
    # Look up the Entrez ID for each symbol and keep only those that were mapped
    unname(symbolToEntrez[as.character(complex)])
  })
  
  # Clean out any NAs that resulted from unmapped symbols
  myComplexesEntrez <- lapply(myComplexesEntrez, function(x) x[!is.na(x)])

  # Inspect the result
  print(myComplexesEntrez)
}
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## $CPLX1
    ## [1] "5430" "5431" "5432"
    ## 
    ## $CPLX2
    ## [1] "983"  "891"  "9133"

Now, the `myComplexesEntrez` list is directly compatible with the
`goGmt` generated from Bioconductor.

### Method 4: From Reactome Pathways via Bioconductor

Similarly, the `reactome.db` package provides pathway annotations, which
also use Entrez IDs.

``` r

if (requireNamespace("reactome.db", quietly = TRUE) && 
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages(library(reactome.db))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  reactomeGmt <- getReactomeGmt(speciesDb = org.Hs.eg.db)

  # Inspect the identifiers
  reactomeGmt[1]
}
```

    ## Fetching Reactome pathway gene sets from Bioconductor...

    ## $`Homo sapiens: Hemostasis`
    ##   [1] "1"      "10019"  "10112"  "10125"  "10125"  "1017"   "10184"  "1020"  
    ##   [9] "10235"  "10235"  "10242"  "10257"  "10335"  "10362"  "10376"  "10381" 
    ##  [17] "10382"  "10383"  "10411"  "10447"  "10451"  "10451"  "10461"  "1048"  
    ##  [25] "10487"  "10490"  "10544"  "10603"  "1062"   "10630"  "10672"  "10672" 
    ##  [33] "10681"  "10681"  "1072"   "10749"  "1084"   "10846"  "1088"   "10916" 
    ##  [41] "10938"  "11004"  "11069"  "11093"  "11127"  "11127"  "11136"  "11216" 
    ##  [49] "112714" "11311"  "113220" "11343"  "113457" "1191"   "12"     "120425"
    ##  [57] "124602" "126961" "1277"   "1277"   "1278"   "1278"   "135228" "139189"
    ##  [65] "139322" "1398"   "139818" "14"     "140628" "140885" "1432"   "1445"  
    ##  [73] "146850" "146909" "147700" "150"    "151"    "152"    "1521"   "1525"  
    ##  [81] "1606"   "1607"   "1608"   "160851" "1609"   "161882" "1675"   "1793"  
    ##  [89] "1794"   "1795"   "1893"   "1950"   "197"    "2"      "207"    "213"   
    ##  [97] "2147"   "2149"   "2151"   "2151"   "2152"   "2153"   "2155"   "2157"  
    ## [105] "2158"   "2159"   "2160"   "2161"   "2162"   "2165"   "2207"   "221037"
    ## [113] "221458" "221955" "2243"   "2244"   "226"    "2266"   "2268"   "2277"  
    ## [121] "22915"  "22920"  "22920"  "22927"  "22953"  "23028"  "23046"  "23052" 
    ## [129] "23095"  "2316"   "23186"  "23303"  "23348"  "2335"   "23414"  "23428" 
    ## [137] "23468"  "23533"  "23539"  "23657"  "23759"  "23764"  "23764"  "24137" 
    ## [145] "25"     "2534"   "259215" "25942"  "25970"  "26090"  "26153"  "2621"  
    ## [153] "2623"   "2624"   "2625"   "2626"   "2627"   "27040"  "27094"  "27154" 
    ## [161] "2734"   "27345"  "2767"   "2767"   "2768"   "2768"   "2769"   "2769"  
    ## [169] "2770"   "2771"   "2773"   "2776"   "2776"   "2778"   "2778"   "2782"  
    ## [177] "2782"   "2783"   "2783"   "2784"   "2784"   "2785"   "2785"   "2786"  
    ## [185] "2786"   "2787"   "2787"   "2788"   "2788"   "2790"   "2790"   "2791"  
    ## [193] "2791"   "2792"   "2792"   "2793"   "2793"   "2811"   "2812"   "2814"  
    ## [201] "2815"   "2817"   "284"    "284"    "285"    "285643" "2885"   "2886"  
    ## [209] "2888"   "29106"  "29127"  "2977"   "29789"  "29802"  "2982"   "2983"  
    ## [217] "2993"   "2994"   "2995"   "302"    "3020"   "3021"   "3043"   "3045"  
    ## [225] "3046"   "3047"   "3048"   "3053"   "3065"   "3066"   "308"    "3082"  
    ## [233] "30845"  "30846"  "3265"   "3273"   "3309"   "333932" "334"    "335"   
    ## [241] "338"    "3439"   "3440"   "3441"   "3442"   "3443"   "3444"   "3445"  
    ## [249] "3446"   "3447"   "3448"   "3449"   "3451"   "3452"   "3456"   "346562"
    ## [257] "347688" "347733" "3479"   "3481"   "350"    "351"    "3512"   "3543"  
    ## [265] "3635"   "3655"   "3659"   "3660"   "3671"   "3672"   "3673"   "3674"  
    ## [273] "3675"   "3676"   "3678"   "3683"   "3684"   "3685"   "3687"   "3688"  
    ## [281] "3688"   "3689"   "3690"   "3699"   "3700"   "3705"   "3708"   "3709"  
    ## [289] "3710"   "3717"   "374354" "3778"   "3779"   "3796"   "3797"   "3797"  
    ## [297] "3798"   "3799"   "3818"   "3827"   "3831"   "3832"   "3833"   "3833"  
    ## [305] "3834"   "3835"   "3845"   "387"    "388"    "3897"   "391"    "3920"  
    ## [313] "3932"   "3937"   "3959"   "4067"   "4072"   "408"    "409"    "4097"  
    ## [321] "4097"   "4099"   "4267"   "4282"   "4312"   "4352"   "440533" "4602"  
    ## [329] "462"    "4680"   "4778"   "4778"   "481"    "482"    "483"    "4842"  
    ## [337] "4843"   "4846"   "487"    "488"    "489"    "4893"   "490"    "491"   
    ## [345] "492"    "493"    "4973"   "5004"   "5005"   "5023"   "5024"   "5025"  
    ## [353] "5026"   "5027"   "5028"   "5028"   "5051"   "5054"   "5055"   "50808" 
    ## [361] "50848"  "50940"  "5099"   "5104"   "51097"  "51206"  "51266"  "51317" 
    ## [369] "5136"   "51368"  "51378"  "5138"   "5152"   "5153"   "5154"   "5155"  
    ## [377] "51571"  "5170"   "5170"   "51706"  "51744"  "5175"   "51764"  "51764" 
    ## [385] "51807"  "5196"   "5197"   "5216"   "5265"   "5267"   "5269"   "5270"  
    ## [393] "5271"   "5290"   "5291"   "5294"   "5295"   "5296"   "5321"   "5327"  
    ## [401] "5328"   "5329"   "5335"   "5336"   "5340"   "5341"   "5345"   "54210" 
    ## [409] "54331"  "54331"  "54495"  "54518"  "54676"  "547"    "5473"   "5478"  
    ## [417] "5478"   "54863"  "55083"  "5515"   "5516"   "5518"   "5519"   "5525"  
    ## [425] "5526"   "5527"   "5528"   "5529"   "55423"  "5547"   "5552"   "55582" 
    ## [433] "55604"  "55605"  "55614"  "55619"  "5566"   "55664"  "55669"  "5567"  
    ## [441] "5568"   "5573"   "5575"   "5576"   "5577"   "5578"   "5579"   "5580"  
    ## [449] "5581"   "5582"   "5583"   "5590"   "5592"   "5592"   "5593"   "5593"  
    ## [457] "5594"   "5595"   "55970"  "55970"  "5624"   "5627"   "56301"  "5657"  
    ## [465] "5660"   "5669"   "5670"   "5671"   "5672"   "5673"   "5675"   "5676"  
    ## [473] "5678"   "5680"   "56992"  "56992"  "57113"  "57126"  "5739"   "5739"  
    ## [481] "57406"  "5747"   "57572"  "5768"   "5770"   "5777"   "5781"   "58494" 
    ## [489] "5868"   "5874"   "5879"   "5880"   "5889"   "5890"   "5894"   "5906"  
    ## [497] "5908"   "5919"   "59345"  "59345"  "60"     "6281"   "634"    "6382"  
    ## [505] "6383"   "6385"   "6401"   "6402"   "6403"   "6404"   "6414"   "64145" 
    ## [513] "64147"  "64147"  "6464"   "64780"  "64805"  "64837"  "6520"   "653604"
    ## [521] "6543"   "6546"   "6547"   "6566"   "66005"  "6647"   "6654"   "6678"  
    ## [529] "6693"   "6694"   "6714"   "6717"   "6786"   "6810"   "6813"   "6814"  
    ## [537] "682"    "6850"   "6850"   "6915"   "7010"   "7010"   "7018"   "7035"  
    ## [545] "7040"   "7042"   "7043"   "7044"   "7056"   "7057"   "7066"   "7076"  
    ## [553] "7078"   "708"    "7094"   "710"    "7102"   "7114"   "7123"   "7157"  
    ## [561] "7222"   "7225"   "7273"   "7277"   "7278"   "7280"   "7409"   "7409"  
    ## [569] "7410"   "7410"   "7414"   "7422"   "7423"   "7424"   "7441"   "7450"  
    ## [577] "7465"   "747"    "7525"   "7534"   "7804"   "7846"   "7873"   "78991" 
    ## [585] "7975"   "7975"   "79861"  "80005"  "801"    "80228"  "805"    "80739" 
    ## [593] "808"    "81"     "81027"  "813"    "813949" "8140"   "8165"   "81704" 
    ## [601] "81928"  "81930"  "81930"  "829"    "829"    "830"    "830"    "832"   
    ## [609] "832"    "8350"   "8351"   "8352"   "8353"   "8354"   "8355"   "8356"  
    ## [617] "8357"   "8358"   "83692"  "83700"  "83706"  "83953"  "8407"   "84617" 
    ## [625] "84643"  "84790"  "84876"  "8503"   "8515"   "8525"   "8526"   "8527"  
    ## [633] "85440"  "857"    "8654"   "87"     "8793"   "8795"   "8797"   "88"    
    ## [641] "8832"   "8968"   "89953"  "9002"   "9002"   "9046"   "9056"   "9057"  
    ## [649] "90952"  "90990"  "9123"   "9127"   "914"    "9162"   "91768"  "928"   
    ## [657] "9371"   "9371"   "94121"  "94235"  "94235"  "9463"   "948"    "9493"  
    ## [665] "9564"   "9585"   "960"    "961"    "962"    "9630"   "9630"   "965"   
    ## [673] "967"    "9672"   "972"    "9732"   "9749"   "9927"   "9948"   "998"

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
to virtually any organism for which you have complex and annotation
data.

Always check your identifiers before running the main analysis functions
to ensure a smooth and successful workflow.
