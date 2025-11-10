#' Example Protein Complex List
#'
#' @description A dataset containing a list of 622 putative protein complexes
#' identified in the study using an integrative global proteomic profiling
#' approach. This dataset is intended for use in package examples and vignettes.
#'
#' @details These complexes were derived from chromatographic separation of
#' cultured human cell extracts into more than one thousand biochemical fractions,
#' followed by quantitative tandem mass spectrometry analysis. Most of the 622 
#' complexes are linked to core biological processes and include candidate 
#' disease genes and unannotated proteins.
#'
#' @format A named list with 622 elements. Each element is a character vector
#' of protein identifiers representing a single complex.
#'
#' @source Supplementary Table 3 (first tab) from:
#' Havugimana PC, Hart GT, Nepusz T, Yang H, Turinsky AL, Li Z, Wang PI, Boutz DR,
#' Fong V, Phanse S, Babu M, Craig SA, Hu P, Wan C, Vlasblom J, Dar VU, Bezginov A,
#' Clark GW, Wu GC, Wodak SJ, Tillier ER, Paccanaro A, Marcotte EM, Emili A.
#' *A census of human soluble protein complexes.* Cell. 2012 Aug 31;150(5):1068-81.
#' doi: 10.1016/j.cell.2012.08.011. PMID: 22939629; PMCID: PMC3477804.
#'
#' @seealso \code{\link{referenceComplexes}}
"demoComplexes"

#' Reference Protein Complex Set (CORUM)
#'
#' @description A reference dataset containing 324 merged CORUM protein complexes
#' used in training protein-protein interaction scoring and clustering optimization
#' procedures in the study.
#'
#' @details These complexes were curated from the CORUM database and merged to
#' reduce redundancy for benchmarking purposes.
#'
#' @format A named list with 324 elements. Each element is a character vector
#' of protein identifiers.
#'
#' @source Supplementary Table 3 (second tab) from:
#' Havugimana PC, Hart GT, Nepusz T, Yang H, Turinsky AL, Li Z, Wang PI, Boutz DR,
#' Fong V, Phanse S, Babu M, Craig SA, Hu P, Wan C, Vlasblom J, Dar VU, Bezginov A,
#' Clark GW, Wu GC, Wodak SJ, Tillier ER, Paccanaro A, Marcotte EM, Emili A.
#' *A census of human soluble protein complexes.* Cell. 2012 Aug 31;150(5):1068-81.
#' doi: 10.1016/j.cell.2012.08.011. PMID: 22939629; PMCID: PMC3477804.
#'
#' @seealso \code{\link{demoComplexes}}
"referenceComplexes"
