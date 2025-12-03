utils::globalVariables(c("primaryFunctionalDomain", "domainSize", "x", "y",
                         "weight", "sizeMapping", "colorHex", "proteins",
                         "source_complex_id", "target_complex_id", "color",
                         "proteinCount", ".data", "score"))

# ============================================================================
# INTERNAL HELPER FUNCTIONS
# ============================================================================

#' Build Base Graph Structure (Safeguarded)
#' @keywords internal
.build_base_graph <- function(layoutDf, edgesDf, edgeColor, verbose) {
  layoutDf <- as.data.frame(layoutDf)
  edgesDf <- as.data.frame(edgesDf)
  
  if (nrow(edgesDf) == 0) {
    graphForPlotting <- tidygraph::tbl_graph(nodes = layoutDf, directed = FALSE)
    p <- ggraph::ggraph(graphForPlotting, layout = 'manual', x = x, y = y)
  } else {
    graphForPlotting <- tidygraph::tbl_graph(nodes = layoutDf, edges = edgesDf, directed = FALSE)
    p <- ggraph::ggraph(graphForPlotting, layout = 'manual', x = x, y = y) +
      ggraph::geom_edge_fan(
        ggplot2::aes(alpha = weight), color = edgeColor, width = 0.25
      )
  }
  return(p)
}

#' Apply Color Mapping (Quantitative or Categorical)
#' @keywords internal
.apply_color_mapping <- function(p, layoutDf, color.by, color.palette, 
                                 color.legend.title, geom_type = c("color", "fill"),
                                 unenrichedColor = "#CCCCCC") {
  geom_type <- match.arg(geom_type)
  
  if (is.null(color.by)) {
    # Categorical mode (Identity or Manual)
    # Handled within the main functions usually, but we can standardize here if needed.
    # For DirectLabels, we use identity. For Legend, we use manual.
    return(p)
  }
  
  # Continuous mode
  if (!color.by %in% names(layoutDf)) {
    stop(sprintf("Column '%s' not found in node data.", color.by))
  }
  
  legend_title <- if (!is.null(color.legend.title)) color.legend.title else color.by
  
  scale_fn <- if (geom_type == "color") {
    list(viridis = ggplot2::scale_color_viridis_c,
         gradientn = ggplot2::scale_color_gradientn)
  } else {
    list(viridis = ggplot2::scale_fill_viridis_c,
         gradientn = ggplot2::scale_fill_gradientn)
  }
  
  if (length(color.palette) == 1 && 
      color.palette %in% c("viridis", "magma", "inferno", "plasma", "cividis")) {
    p <- p + scale_fn$viridis(option = color.palette, name = legend_title, na.value = "grey50")
  } else {
    p <- p + scale_fn$gradientn(colors = color.palette, name = legend_title, na.value = "grey50")
  }
  
  return(p)
}

#' Create Tooltip HTML for visNetwork
#' @keywords internal
.create_tooltip <- function(complexId, primaryFunctionalDomain, proteinCount, 
                            proteins, score = NULL, color.by = NULL, color.value = NULL) {
  base_html <- paste0(
    "<div style='font-family:sans-serif; text-align:left;'>",
    "<b>Complex:</b> ", complexId, "<br>",
    "<b>Function:</b> ", primaryFunctionalDomain, "<br>"
  )
  
  if (!is.null(score) && !is.na(score)) {
    base_html <- paste0(base_html, "<b>Specificity Score:</b> ", round(score, 2), "<br>")
  }
  
  # MODIFIED: Check for NA to handle placeholder values correctly.
  if (!is.null(color.by) && !is.na(color.by) && !is.null(color.value) && !is.na(color.value)) {
    base_html <- paste0(
      base_html,
      "<b>", color.by, ":</b> ", round(color.value, 2), "<br>"
    )
  }
  
  paste0(
    base_html,
    "<b>Protein Count:</b> ", proteinCount, "<hr>",
    "<b>Members:</b><br>", gsub(",", ", ", proteins), "</div>"
  )
}


#' Create Label Data for Direct Labels
#' @keywords internal
.create_label_data <- function(layoutDf, centroid_threshold = 2) {
  if (!"primaryFunctionalDomain" %in% names(layoutDf)) return(list(centroid=data.frame(), direct=data.frame()))
  
  labelData <- layoutDf %>%
    dplyr::filter(primaryFunctionalDomain != "Unenriched" & !is.na(primaryFunctionalDomain)) %>%
    dplyr::add_count(primaryFunctionalDomain, name = "domainSize")
  
  if (nrow(labelData) == 0) return(list(centroid = data.frame(), direct = data.frame()))
  
  centroidLabels <- labelData %>%
    dplyr::filter(domainSize > centroid_threshold) %>%
    dplyr::group_by(primaryFunctionalDomain) %>%
    dplyr::summarise(x = mean(x, na.rm=TRUE), y = mean(y, na.rm=TRUE), .groups = "drop")
  
  directLabels <- labelData %>% 
    dplyr::filter(domainSize <= centroid_threshold)
  
  list(centroid = as.data.frame(centroidLabels), direct = as.data.frame(directLabels))
}

# ============================================================================
# EXPORTED VISUALIZATION FUNCTIONS
# ============================================================================

#' Visualize a Complex Map with Direct Node Labels
#'
#' @description
#' Creates a static visualization where functional labels are placed directly on the plot.
#'
#' @param layoutDf A data frame containing node attributes and layout coordinates.
#' @param edgesDf A data frame containing the network edges.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param bgColor Background color.
#' @param edgeColor Edge color.
#' @param nodeSizeRange Min and max node size.
#' @param labelFillColor Background fill color for labels.
#' @param centroid_threshold Domains with > N complexes get a single centroid label.
#' @param color.by Optional numeric column for continuous coloring.
#' @param color.palette Palette name ("viridis", "plasma") or vector of colors.
#' @param color.legend.title Title for color legend.
#' @param verbose Logical.
#'
#' @return A `ggplot` object.
#'
#' @export
visualizeMapDirectLabels <- function(
    layoutDf, edgesDf,
    title = "ComplexMap Functional Landscape",
    subtitle = "Nodes colored by Specificity-Weighted Function",
    bgColor = "black", edgeColor = "white",
    nodeSizeRange = c(2, 10), 
    labelFillColor = ggplot2::alpha("white", 0.7),
    centroid_threshold = 2,
    color.by = NULL,
    color.palette = "viridis",
    color.legend.title = NULL,
    verbose = TRUE) {
  
  if (verbose) message("Visualizing diverse landscape with direct labels...")
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.", call. = FALSE)
  }
  
  p <- .build_base_graph(layoutDf, edgesDf, edgeColor, verbose)
  
  # Logic: If color.by is NULL, use calculated 'colorHex'. If present, use that column.
  if (is.null(color.by)) {
    p <- p + 
      ggraph::geom_node_point(
        ggplot2::aes(size = sizeMapping, color = colorHex)
      ) +
      ggplot2::scale_color_identity()
  } else {
    p <- p + 
      ggraph::geom_node_point(
        ggplot2::aes(size = sizeMapping, color = .data[[color.by]])
      )
    p <- .apply_color_mapping(p, layoutDf, color.by, color.palette, color.legend.title, geom_type="color")
  }
  
  # Safe Label Addition
  labels <- .create_label_data(layoutDf, centroid_threshold)
  
  if (nrow(labels$centroid) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labels$centroid,
      ggplot2::aes(x = x, y = y, label = stringr::str_wrap(primaryFunctionalDomain, 20)),
      fill = labelFillColor, fontface = 'bold', color = "black",
      segment.color = edgeColor
    )
  }
  
  if (nrow(labels$direct) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = labels$direct,
      ggplot2::aes(x = x, y = y, label = stringr::str_wrap(primaryFunctionalDomain, 20)),
      fill = labelFillColor, fontface = 'bold', color = "black",
      segment.color = edgeColor, box.padding = 0.5
    )
  }
  
  p +
    ggplot2::scale_size_continuous(
      range = nodeSizeRange, name = "Size (log2)",
      guide = ggplot2::guide_legend(override.aes = list(color = "white"))
    ) +
    ggraph::scale_edge_alpha_continuous(range = c(0.1, 1), name = "Similarity") +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bgColor),
      panel.background = ggplot2::element_rect(fill = bgColor),
      plot.title = ggplot2::element_text(color = edgeColor, hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(color = edgeColor, hjust = 0.5),
      legend.position = "right",
      legend.text = ggplot2::element_text(color = edgeColor),
      legend.title = ggplot2::element_text(color = edgeColor),
      legend.background = ggplot2::element_rect(fill = bgColor)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
}

#' Visualize a Complex Map with a Color Legend
#'
#' @description
#' Creates a static visualization with a discrete legend for functional domains.
#'
#' @inheritParams visualizeMapDirectLabels
#' @param unenrichedColor Color for unenriched nodes.
#'
#' @return A `ggplot` object.
#'
#' @export
visualizeMapWithLegend <- function(
    layoutDf, edgesDf,
    title = "ComplexMap Functional Landscape",
    subtitle = "Nodes colored by Specificity-Weighted Function",
    bgColor = "black", edgeColor = "white",
    nodeSizeRange = c(2, 10), 
    unenrichedColor = "#CCCCCC",
    color.by = NULL,
    color.palette = "plasma",
    color.legend.title = NULL,
    verbose = TRUE) {
  
  if (verbose) message("Visualizing diverse landscape with color legend...")
  
  p <- .build_base_graph(layoutDf, edgesDf, edgeColor, verbose)
  
  if (is.null(color.by)) {
    # Categorical Functional Domains
    if (nrow(layoutDf) > 0 && "primaryFunctionalDomain" %in% names(layoutDf)) {
      domainColorMap <- layoutDf %>%
        dplyr::filter(primaryFunctionalDomain != "Unenriched" & !is.na(primaryFunctionalDomain))
      
      if (nrow(domainColorMap) > 0) {
        domainColorMap <- domainColorMap %>%
          dplyr::group_by(primaryFunctionalDomain) %>%
          dplyr::summarise(color = names(which.max(table(colorHex))), .groups = "drop")
        legendColors <- tibble::deframe(domainColorMap)
      } else {
        legendColors <- character(0)
      }
      legendColors["Unenriched"] <- unenrichedColor
      
      p <- p + 
        ggraph::geom_node_point(
          ggplot2::aes(fill = primaryFunctionalDomain, size = sizeMapping),
          shape = 21, color = edgeColor, stroke = 0.2
        ) +
        ggplot2::scale_fill_manual(
          name = "Functional Domain", values = legendColors,
          guide = ggplot2::guide_legend(override.aes = list(size = 5), ncol = 2)
        )
    } else {
      p <- p + ggraph::geom_node_point(ggplot2::aes(size = sizeMapping), color = unenrichedColor)
    }
  } else {
    # Continuous Coloring
    p <- p + 
      ggraph::geom_node_point(
        ggplot2::aes(fill = .data[[color.by]], size = sizeMapping),
        shape = 21, color = edgeColor, stroke = 0.2
      )
    p <- .apply_color_mapping(p, layoutDf, color.by, color.palette, color.legend.title, geom_type="fill")
  }
  
  p +
    ggplot2::scale_size_continuous(range = nodeSizeRange, name = "Size (log2)") +
    ggraph::scale_edge_alpha_continuous(range = c(0.1, 1), name = "Similarity") +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bgColor),
      panel.background = ggplot2::element_rect(fill = bgColor),
      legend.background = ggplot2::element_rect(fill = bgColor),
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(color = edgeColor, size = 8),
      legend.title = ggplot2::element_text(color = edgeColor),
      plot.title = ggplot2::element_text(color = edgeColor, hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(color = edgeColor, hjust = 0.5)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
}

#' Visualize a Complex Map Interactively
#'
#' @description
#' Creates a dynamic `visNetwork` widget.
#'
#' @inheritParams visualizeMapDirectLabels
#' @param width Widget width.
#' @param height Widget height.
#' @param physicsEnabled Enable physics simulation.
#'
#' @return A `visNetwork` widget.
#'
#' @export
visualizeMapInteractive <- function(
    layoutDf, edgesDf,
    width = "100%", height = "90vh",
    title = "ComplexMap Functional Landscape",
    physicsEnabled = FALSE,
    color.by = NULL,
    color.palette = "viridis",
    verbose = TRUE) {
  
  if (verbose) message("Generating interactive visNetwork plot...")
  
  nodesVn <- layoutDf %>%
    dplyr::mutate(
      id = complexId,
      label = complexId,
      value = sizeMapping, 
      x = scales::rescale(x, to = c(-800, 800)),
      y = scales::rescale(-y, to = c(-800, 800))
    )
  
  # Handle color
  if (is.null(color.by)) {
    nodesVn$color <- nodesVn$colorHex
  } else {
    if (!color.by %in% names(layoutDf)) stop(sprintf("Column '%s' not found.", color.by))
    quant_values <- layoutDf[[color.by]]
    # Use scales to map color
    color_mapper <- scales::col_numeric(palette = color.palette, domain = range(quant_values, na.rm=TRUE))
    nodesVn$color <- ifelse(is.na(quant_values), "#888888", color_mapper(quant_values))
  }
  
  # MODIFIED: Prepare arguments for mapply to avoid passing NULL, which causes
  # the function to return a zero-length list and breaks the assignment.
  score_arg <- if ("score" %in% names(nodesVn)) nodesVn$score else NA
  color_by_name_arg <- if (!is.null(color.by)) color.by else NA
  color_value_arg <- if (!is.null(color.by) && color.by %in% names(nodesVn)) nodesVn[[color.by]] else NA
  
  nodesVn$title <- mapply(
    FUN = .create_tooltip,
    nodesVn$complexId,
    nodesVn$primaryFunctionalDomain,
    nodesVn$proteinCount,
    nodesVn$proteins,
    score = score_arg,
    color.by = color_by_name_arg,
    color.value = color_value_arg,
    SIMPLIFY = FALSE
  )
  
  edgesVn <- edgesDf %>%
    dplyr::mutate(from = source_complex_id, to = target_complex_id, value = weight)
  
  visNetwork::visNetwork(nodes = nodesVn, edges = edgesVn, width = width, height = height, main = title) %>%
    visNetwork::visNodes(shape = "dot", borderWidth = 2, shadow = TRUE) %>%
    visNetwork::visEdges(smooth = FALSE, color = list(color = "#888888", highlight = "#00BFFF")) %>%
    visNetwork::visPhysics(enabled = physicsEnabled) %>%
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE,
      selectedBy = "primaryFunctionalDomain"
    ) %>%
    visNetwork::visLayout(randomSeed = 123) %>%
    visNetwork::visInteraction(navigationButtons = TRUE, dragNodes = TRUE, zoomView = TRUE)
}