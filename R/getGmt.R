utils::globalVariables(c("ONTOLOGY", "ENTREZID", "TERM",
                         "gs_name", "gene_symbol"))
#' Fetch Gene Sets from MSigDB
#'
#' @description A wrapper for `msigdbr::msigdbr` to fetch gene sets and format
#' them as a named list (GMT format).
#'
#' @param species The scientific name for the species (e.g., "Homo sapiens").
#' @param collection The MSigDB collection code (e.g., "H" for hallmark).
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A named list where names are gene set names and values are character
#'   vectors of gene symbols.
#'
#' @export
#'
getMsigdbGmt <- function(species="Homo sapiens", collection="H", verbose=TRUE) {
  if (verbose) {
    message("Fetching MSigDB sets (Species: ", species,
            ", Cat: ", collection, ")")
  }
  msigdbr::msigdbr(species=species, collection=collection) %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(genes=list(gene_symbol), .groups='drop') %>%
    tibble::deframe()
}

#' Read a GMT File from a Local Path
#'
#' @description Parses a local GMT (Gene Matrix Transposed) file into the
#' standard named-list format.
#'
#' @param filepath The path to the local .gmt file.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A named list where names are gene set names and values are character
#'   vectors of genes.
#'
#' @export
#'
getGmtFromFile <- function(filepath, verbose=TRUE) {
  if (verbose) message("Fetching gene sets from local file: ", filepath)
  gmtLines <- readLines(filepath)
  gmtList <- lapply(gmtLines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    genes <- parts[-c(1, 2)]
    return(genes[genes != ""])
  })
  names(gmtList) <- vapply(gmtLines, function(line) {
    strsplit(line, "\t")[[1]][1]
  }, character(1))
  return(gmtList)
}

#' Fetch Gene Ontology (GO) Gene Sets
#'
#' @description Fetches GO terms and their associated genes from a Bioconductor
#' AnnotationDb package.
#'
#' @param speciesDb An AnnotationDb object for the target species (e.g.,
#'   `org.Hs.eg.db`).
#' @param ontology The GO ontology to fetch. One of "BP", "MF", or "CC".
#' @param minGmtSize The minimum number of genes for a GO term to be included.
#' @param maxGmtSize The maximum number of genes for a GO term to be included.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A named list where names are GO terms and values are character
#'   vectors of Entrez IDs.
#'
#' @export
#'
getGoGmt <- function(speciesDb, ontology="BP", minGmtSize=10,
                     maxGmtSize=500, verbose=TRUE) {
  if (verbose) {
    message("Fetching GO (", ontology, ") gene sets from Bioconductor...")
  }
  goKeys <- AnnotationDbi::keys(speciesDb, keytype="GO")
  go_to_entrez <- AnnotationDbi::select(
    speciesDb, keys=goKeys, columns=c("ENTREZID", "ONTOLOGY"), keytype="GO"
  ) %>%
    dplyr::filter(ONTOLOGY == ontology)
  
  go_to_term <- AnnotationDbi::select(
    GO.db::GO.db, keys=unique(go_to_entrez$GO),
    columns="TERM", keytype="GOID"
  )
  
  gmtUnfiltered <- go_to_entrez %>%
    dplyr::left_join(go_to_term, by=c("GO"="GOID")) %>%
    dplyr::filter(!is.na(ENTREZID) & !is.na(TERM)) %>%
    dplyr::group_by(TERM) %>%
    dplyr::summarise(genes=list(unique(ENTREZID)), .groups='drop') %>%
    tibble::deframe()
  
  gmtFiltered <- Filter(function(genes) {
    num_genes <- length(genes)
    num_genes >= minGmtSize && num_genes <= maxGmtSize
  }, gmtUnfiltered)
  
  if (verbose) {
    message("Retained ", length(gmtFiltered), " GO terms (size ",
            minGmtSize, "-", maxGmtSize, ") out of ", length(gmtUnfiltered))
  }
  return(gmtFiltered)
}

#' Fetch Reactome Pathway Gene Sets
#'
#' @description Fetches Reactome pathways and their associated genes from
#' the `reactome.db` Bioconductor annotation package.
#'
#' @param speciesDb An AnnotationDb object for the target species (e.g.,
#'   `org.Hs.eg.db`). This is used to filter pathways for the correct species.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A named list where names are Reactome pathway names and values are
#'   character vectors of Entrez IDs.
#'
#' @export
#'
getReactomeGmt <- function(speciesDb, verbose=TRUE) {
  if (verbose) {
    message("Fetching Reactome pathway gene sets from Bioconductor...")
  }
  reactomeToEntrez <- as.list(reactome.db::reactomePATHID2EXTID)
  reactomeToName <- as.list(reactome.db::reactomePATHID2NAME)
  
  meta <- AnnotationDbi::metadata(speciesDb)
  speciesName <- meta[meta$name == "ORGANISM", "value"]
  
  gmt <- lapply(names(reactomeToEntrez), function(pathId) {
    pathwayName <- reactomeToName[[pathId]]
    isCorrectSpecies <- !is.null(pathwayName) &&
      grepl(speciesName, pathwayName, ignore.case=TRUE)
    
    if (isCorrectSpecies) {
      return(reactomeToEntrez[[pathId]])
    }
    return(NULL)
  })
  
  names(gmt) <- vapply(names(reactomeToEntrez), function(id) {
    name <- reactomeToName[[id]]
    if (is.null(name)) "" else name
  }, character(1))
  
  # Remove the NULL entries for pathways from other species
  Filter(Negate(is.null), gmt)
}

#' Get the Path to an Example GMT File
#'
#' @description Provides the full system path to an example gene set (GMT)
#' file included with the ComplexMap package.
#'
#' @details The included file is the BioCarta gene set collection from the
#' Molecular Signatures Database (MSigDB v2025.1). This function makes it easy
#' to access the file for use in examples and vignettes.
#'
#' @return A character string containing the full path to the example GMT file.
#'
#' @export
#' @examples
#' # Get the path
#' gmtPath <- getExampleGmt()
#'
#' # You can then read the file using the path
#' if (file.exists(gmtPath)) {
#'   exampleGmt <- getGmtFromFile(gmtPath)
#' }
#'
getExampleGmt <- function() {
  system.file("extdata", "c2.cp.biocarta.v2025.1.Hs.symbols.gmt",
              package="ComplexMap", mustWork=TRUE)
}