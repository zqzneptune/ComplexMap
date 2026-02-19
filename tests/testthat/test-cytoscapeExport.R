## tests/testthat/test-cytoscapeExport.R
library(testthat)
library(ComplexMap)
library(jsonlite)

# ── Minimal mock data setup ──────────────────────────────────────────────────

make_mock_cmap <- function() {
  nodes <- tibble::tibble(
    complexId               = c("Cpx1", "Cpx2", "Cpx3"),
    proteinCount            = c(5L, 8L, 3L),
    proteins                = c("A,B,C,D,E", "B,C,D,E,F,G,H,I", "G,H,I"),
    primaryFunctionalDomain = c("Cell Cycle", "DNA Repair", "Unenriched"),
    topEnrichedFunctions    = c("G1/S transition", "Nucleotide excision repair", NA_character_),
    colorHex                = c("#1C7ED6", "#F03E3E", "#CCCCCC"),
    x                       = c(-1.5,  0.3,  2.1),
    y                       = c( 0.8, -1.2,  0.1),
    betweenness             = c(0.55, 0.80, 0.10),
    degree                  = c(2L, 3L, 1L)
  )

  edges <- tibble::tibble(
    source_complex_id = c("Cpx1", "Cpx2"),
    target_complex_id = c("Cpx2", "Cpx3"),
    compSim           = c(0.42, 0.28),
    funcSim           = c(0.35, 0.10),
    weight            = c(0.40, 0.25),
    similarity_mode   = c("combined", "combined")
  )

  structure(
    list(
      nodes = nodes,
      edges = edges,
      layout_info = list(method = "fr", seed = 123, scaling = "raw")
    ),
    class = "ComplexMap"
  )
}

# ── toCytoscapeJSON ──────────────────────────────────────────────────────────

test_that("toCytoscapeJSON returns a valid JSON string", {
  cmap <- make_mock_cmap()
  json_str <- toCytoscapeJSON(cmap, prune = FALSE, verbose = FALSE)

  expect_type(json_str, "character")
  expect_true(nchar(json_str) > 10)

  parsed <- jsonlite::fromJSON(json_str, simplifyVector = FALSE)
  expect_type(parsed, "list")
})

test_that("node count in JSON matches ComplexMap node table", {
  cmap <- make_mock_cmap()
  json_str <- toCytoscapeJSON(cmap, prune = FALSE, verbose = FALSE)
  parsed <- jsonlite::fromJSON(json_str, simplifyVector = FALSE)

  # JSON has both nodes (no source) and edges (have source)
  nodes_in_json <- Filter(function(el) is.null(el$data$source), parsed)
  expect_equal(length(nodes_in_json), nrow(getNodeTable(cmap)))
})

test_that("edge count in JSON matches ComplexMap edge table (when not pruned)", {
  cmap <- make_mock_cmap()
  json_str <- toCytoscapeJSON(cmap, prune = FALSE, verbose = FALSE)
  parsed <- jsonlite::fromJSON(json_str, simplifyVector = FALSE)

  edges_in_json <- Filter(function(el) !is.null(el$data$source), parsed)
  expect_equal(length(edges_in_json), nrow(getEdgeTable(cmap)))
})

test_that("node data contains required Cytoscape fields", {
  cmap <- make_mock_cmap()
  json_str <- toCytoscapeJSON(cmap, prune = FALSE, verbose = FALSE)
  parsed <- jsonlite::fromJSON(json_str, simplifyVector = FALSE)

  node_el <- Filter(function(el) is.null(el$data$source), parsed)[[1]]
  required <- c("id", "label", "proteinCount", "primaryFunctionalDomain",
                "colorHex", "betweenness", "degree")
  for (field in required) {
    expect_true(!is.null(node_el$data[[field]]),
                label = paste("Missing field:", field))
  }
})

test_that("node positions are within [-500, 500] after scaling", {
  cmap <- make_mock_cmap()
  json_str <- toCytoscapeJSON(cmap, prune = FALSE, verbose = FALSE)
  parsed <- jsonlite::fromJSON(json_str, simplifyVector = FALSE)

  nodes_in_json <- Filter(function(el) is.null(el$data$source), parsed)
  xs <- sapply(nodes_in_json, function(el) el$position$x)
  ys <- sapply(nodes_in_json, function(el) el$position$y)

  expect_true(all(xs >= -501) && all(xs <= 501))
  expect_true(all(ys >= -501) && all(ys <= 501))
})

# ── scaleLayoutForCytoscape ──────────────────────────────────────────────────

test_that("scaleLayoutForCytoscape produces coords within specified range", {
  set.seed(1)
  x <- rnorm(50, mean = 10, sd = 5)
  y <- rnorm(50, mean = -3, sd = 8)

  scaled <- scaleLayoutForCytoscape(x, y, range = c(-500, 500))

  expect_true(all(scaled$x >= -500.1) && all(scaled$x <= 500.1))
  expect_true(all(scaled$y >= -500.1) && all(scaled$y <= 500.1))
})

test_that("scaleLayoutForCytoscape handles single-value input", {
  result <- scaleLayoutForCytoscape(c(5), c(5))
  expect_equal(result$x, 0)
  expect_equal(result$y, 0)
})

# ── pruneNetwork ─────────────────────────────────────────────────────────────

test_that("pruneNetwork top_k reduces edges", {
  cmap <- make_mock_cmap()
  # Both edges should be kept (k=1 still allows 1 per node)
  # but with k=1 on a 3-node/2-edge graph, all edges may survive
  pruned <- pruneNetwork(cmap, method = "top_k", k = 1, verbose = FALSE)
  expect_lte(nrow(getEdgeTable(pruned)), nrow(getEdgeTable(cmap)))
})

test_that("pruneNetwork quantile reduces edges", {
  # Edge weights: 0.40, 0.25 -> quantile(0.75) = 0.3875
  # Only 0.40 survives
  cmap <- make_mock_cmap()
  pruned <- pruneNetwork(cmap, method = "quantile", weight_quantile = 0.75, verbose = FALSE)
  expect_lt(nrow(getEdgeTable(pruned)), nrow(getEdgeTable(cmap)))
})

test_that("pruneNetwork errors on invalid method", {
  cmap <- make_mock_cmap()
  expect_error(pruneNetwork(cmap, method = "invalid", verbose = FALSE))
})

# ── Layout reproducibility ───────────────────────────────────────────────────

test_that("scaleLayoutForCytoscape is deterministic", {
  x <- c(1.1, 2.2, 3.3)
  y <- c(-0.5, 0.0, 0.5)

  r1 <- scaleLayoutForCytoscape(x, y)
  r2 <- scaleLayoutForCytoscape(x, y)

  expect_equal(r1$x, r2$x)
  expect_equal(r1$y, r2$y)
})
