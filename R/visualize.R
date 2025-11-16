utils::globalVariables(c("primaryFunctionalDomain", "domainSize", "x", "y",
                         "weight", "sizeMapping", "colorHex", "proteins",
                         "source_complex_id", "target_complex_id", "color",
                         "proteinCount", ".data"))

# ============================================================================
# INTERNAL HELPER FUNCTIONS
# ============================================================================

#' Build Base Graph Structure
#' @keywords internal
.build_base_graph <- function(layoutDf, edgesDf, edgeColor, verbose) {
  if (verbose) message("Building base graph structure...")
  
  graphForPlotting <- tidygraph::tbl_graph(
    nodes = layoutDf, edges = edgesDf, directed = FALSE
  )
  
  ggraph::ggraph(graphForPlotting, layout = 'manual', x = x, y = y) +
    ggraph::geom_edge_fan(
      ggplot2::aes(alpha = weight), color = edgeColor, width = 0.25
    )
}

#' Apply Color Mapping (Categorical or Continuous)
#' @keywords internal
.apply_color_mapping <- function(p, layoutDf, color.by, color.palette, 
                                 color.legend.title, geom_type = c("color", "fill"),
                                 unenrichedColor = "#CCCCCC") {
  geom_type <- match.arg(geom_type)
  
  if (is.null(color.by)) {
    # Categorical mode
    if (geom_type == "color") {
      p <- p + ggplot2::scale_color_identity()
    } else {
      # Build domain color map
      domainColorMap <- layoutDf %>%
        dplyr::filter(
          primaryFunctionalDomain != "Unenriched" &
            !is.na(primaryFunctionalDomain)
        ) %>%
        dplyr::group_by(primaryFunctionalDomain) %>%
        dplyr::summarise(
          color = names(which.max(table(colorHex))), .groups = "drop"
        )
      
      legendColors <- tibble::deframe(domainColorMap)
      legendColors["Unenriched"] <- unenrichedColor
      
      p <- p + 
        ggplot2::scale_fill_manual(name = "Functional Domain", values = legendColors) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5)))
    }
  } else {
    # Continuous mode
    if (!color.by %in% names(layoutDf)) {
      stop(sprintf("Column '%s' not found in the provided node data frame.", color.by))
    }
    
    legend_title <- 
      if (!is.null(color.legend.title)) color.legend.title else color.by
    
    scale_fn <- if (geom_type == "color") {
      list(viridis = ggplot2::scale_color_viridis_c,
           gradientn = ggplot2::scale_color_gradientn)
    } else {
      list(viridis = ggplot2::scale_fill_viridis_c,
           gradientn = ggplot2::scale_fill_gradientn)
    }
    
    if (length(color.palette) == 1 && 
        color.palette %in% c("viridis", "magma", "inferno", "plasma", "cividis")) {
      p <- p + scale_fn$viridis(
        option = color.palette, 
        name = legend_title, 
        na.value = "grey50"
      )
    } else {
      p <- p + scale_fn$gradientn(
        colors = color.palette, 
        name = legend_title, 
        na.value = "grey50"
      )
    }
  }
  
  return(p)
}

#' Create Label Data for Direct Labels
#' @keywords internal
.create_label_data <- function(layoutDf, centroid_threshold = 2) {
  labelData <- layoutDf %>%
    dplyr::filter(
      primaryFunctionalDomain != "Unenriched" &
        !is.na(primaryFunctionalDomain)
    ) %>%
    dplyr::add_count(primaryFunctionalDomain, name = "domainSize")
  
  centroidLabels <- labelData %>%
    dplyr::filter(domainSize > centroid_threshold) %>%
    dplyr::group_by(primaryFunctionalDomain) %>%
    dplyr::summarise(x = mean(x), y = mean(y), .groups = "drop")
  
  directLabels <- labelData %>% 
    dplyr::filter(domainSize <= centroid_threshold)
  
  list(centroid = centroidLabels, direct = directLabels)
}

#' Apply Common Theme to ggraph Plot
#' @keywords internal
.apply_theme <- function(p, bgColor, edgeColor, fontFamily, 
                         nodeSizeRange, size.legend.title, 
                         title, subtitle) {
  p +
    ggplot2::scale_size_continuous(
      range = nodeSizeRange, 
      name = size.legend.title,
      guide = ggplot2::guide_legend(override.aes = list(color = "white"))
    ) +
    ggraph::scale_edge_alpha_continuous(
      range = c(0.1, 1), name = "Similarity"
    ) +
    ggraph::theme_graph(base_family = fontFamily) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bgColor),
      panel.background = ggplot2::element_rect(fill = bgColor),
      legend.background = ggplot2::element_rect(fill = bgColor),
      legend.key = ggplot2::element_rect(fill = bgColor),
      legend.text = ggplot2::element_text(color = edgeColor),
      legend.title = ggplot2::element_text(color = edgeColor),
      plot.title = ggplot2::element_text(color = edgeColor, hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(color = edgeColor, hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
}

#' Create Tooltip HTML for visNetwork
#' @keywords internal
.create_tooltip <- function(complexId, primaryFunctionalDomain, proteinCount, 
                            proteins, color.by = NULL, color.value = NULL) {
  base_html <- paste0(
    "<div style='font-family:sans-serif; text-align:left;'>",
    "<b>Complex:</b> ", complexId, "<br>",
    "<b>Function:</b> ", primaryFunctionalDomain, "<br>"
  )
  
  if (!is.null(color.by) && !is.null(color.value)) {
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

# ============================================================================
# EXPORTED VISUALIZATION FUNCTIONS
# ============================================================================

#' Visualize a Complex Map with Direct Node Labels
#'
#' @description
#' Creates a static visualization of the complex network using `ggraph`, where
#' functional domain labels are placed directly on the plot near the nodes.
#'
#' @details
#' This function supports two coloring modes:
#' 
#' 1.  **Categorical (default):** When `color.by = NULL`, nodes are colored
#'     using the `colorHex` column, which typically represents functional domains.
#'     No color legend is drawn.
#'     
#' 2.  **Continuous:** When `color.by` is set to the name of a numeric column
#'     in `layoutDf` (e.g., "abundance"), nodes are colored along a continuous
#'     gradient based on that column's values. A color bar legend is drawn.
#'
#' The function places labels at the centroid for large domains (>centroid_threshold 
#' complexes) and directly on nodes for smaller domains.
#'
#' @param layoutDf A data frame containing node attributes and layout
#'   coordinates, typically from `computeMapTopology`.
#' @param edgesDf A data frame containing the network edges.
#' @param title The main title for the plot.
#' @param subtitle The subtitle for the plot.
#' @param bgColor The background color of the plot.
#' @param edgeColor The color of the network edges.
#' @param nodeSizeRange A numeric vector of length 2 specifying the min and max
#'   node size.
#' @param labelFillColor The background fill color for the labels.
#' @param fontFamily The base font family for all plot text.
#' @param size.legend.title The title for the node size legend.
#' @param color.by A character string specifying the name of a numeric column in
#'   `layoutDf` to use for continuous node coloring. If `NULL` (default), the
#'   categorical `colorHex` column is used.
#' @param color.palette A character string or vector of colors for the continuous
#'   gradient (e.g., "viridis", "plasma", or `c("blue", "white", "red")`).
#'   Only used when `color.by` is specified.
#' @param color.legend.title A character string for the title of the continuous
#'   color legend. Defaults to the value of `color.by`.
#' @param centroid_threshold Integer. Domains with more than this many complexes
#'   will use centroid labels. Defaults to 2.
#' @param label_wrap_width Integer. Maximum width for label text wrapping. 
#'   Defaults to 20.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A `ggplot` object representing the network visualization.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # --- Sample Data ---
#' nodes <- tibble::tibble(
#'   complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
#'   primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Metabolism"),
#'   sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#0000FF"),
#'   abundance = c(1.2, -0.5, 0.8)
#' )
#' edges <- tibble::tibble(
#'   source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
#' )
#'
#' # --- Usage 1: Default categorical coloring ---
#' if (requireNamespace("ggrepel", quietly = TRUE)) {
#'   visualizeMapDirectLabels(nodes, edges)
#' }
#'
#' # --- Usage 2: Continuous coloring by 'abundance' ---
#' if (requireNamespace("ggrepel", quietly = TRUE)) {
#'   visualizeMapDirectLabels(nodes, edges, color.by = "abundance")
#' }
#'
visualizeMapDirectLabels <- function(
    layoutDf, edgesDf,
    title = "ComplexMap Functional Landscape",
    subtitle = "Nodes are protein complexes, colored by function",
    bgColor = "black", edgeColor = "white",
    nodeSizeRange = c(2, 10), 
    labelFillColor = ggplot2::alpha("white", 0.7),
    fontFamily = "sans", 
    size.legend.title = "Complex Size (log2)",
    color.by = NULL, 
    color.palette = "viridis", 
    color.legend.title = NULL,
    centroid_threshold = 2,
    label_wrap_width = 20,
    verbose = TRUE) {
  
  if (verbose) message("Visualizing ComplexMap with direct labels...")
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for direct labeling.", call. = FALSE)
  }
  
  # Build base graph
  p <- .build_base_graph(layoutDf, edgesDf, edgeColor, verbose = FALSE)
  
  # Add node layer with appropriate color mapping
  if (is.null(color.by)) {
    p <- p + ggraph::geom_node_point(
      ggplot2::aes(size = sizeMapping, color = colorHex)
    )
  } else {
    p <- p + ggraph::geom_node_point(
      ggplot2::aes(size = sizeMapping, color = .data[[color.by]])
    )
  }
  
  # Apply color scale
  p <- .apply_color_mapping(p, layoutDf, color.by, color.palette, 
                            color.legend.title, geom_type = "color")
  
  # Create and add labels
  labels <- .create_label_data(layoutDf, centroid_threshold)
  
  p <- p +
    ggrepel::geom_label_repel(
      data = labels$centroid,
      ggplot2::aes(x = x, y = y,
                   label = stringr::str_wrap(primaryFunctionalDomain, label_wrap_width)),
      fill = labelFillColor, fontface = 'bold', color = "black",
      segment.color = edgeColor
    ) +
    ggrepel::geom_label_repel(
      data = labels$direct,
      ggplot2::aes(x = x, y = y,
                   label = stringr::str_wrap(primaryFunctionalDomain, label_wrap_width)),
      fill = labelFillColor, fontface = 'bold', color = "black",
      segment.color = edgeColor, box.padding = 0.5, point.padding = 0.5
    )
  
  # Apply theme and return
  .apply_theme(p, bgColor, edgeColor, fontFamily, nodeSizeRange, 
               size.legend.title, title, subtitle)
}


#' Visualize a Complex Map with a Color Legend
#'
#' @description
#' Creates a static visualization of the complex network using `ggraph`, with
#' a legend representing node colors.
#'
#' @details
#' This function supports two coloring modes:
#' 
#' 1.  **Categorical (default):** When `color.by = NULL`, nodes are colored
#'     by their `primaryFunctionalDomain`. A discrete color legend is shown.
#'     
#' 2.  **Continuous:** When `color.by` is set to the name of a numeric column
#'     in `layoutDf`, nodes are colored along a continuous gradient. A color
#'     bar legend is shown.
#'
#' @inheritParams visualizeMapDirectLabels
#' @param unenrichedColor The color for nodes in the "Unenriched" category in
#'   categorical mode.
#'
#' @return A `ggplot` object representing the network visualization.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # --- Sample Data ---
#' nodes <- tibble::tibble(
#'   complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
#'   primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Unenriched"),
#'   sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#CCCCCC"),
#'   purity = c(0.95, 0.87, 0.91)
#' )
#' edges <- tibble::tibble(
#'   source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
#' )
#'
#' # --- Usage 1: Default categorical coloring ---
#' visualizeMapWithLegend(nodes, edges)
#'
#' # --- Usage 2: Continuous coloring ---
#' visualizeMapWithLegend(nodes, edges, color.by = "purity")
#'
visualizeMapWithLegend <- function(
    layoutDf, edgesDf,
    title = "ComplexMap Functional Landscape",
    subtitle = "Nodes are protein complexes, colored by function",
    bgColor = "black", edgeColor = "white",
    nodeSizeRange = c(2, 10), 
    unenrichedColor = "#CCCCCC",
    fontFamily = "sans", 
    size.legend.title = "Complex Size (log2)",
    color.by = NULL, 
    color.palette = "viridis", 
    color.legend.title = NULL,
    verbose = TRUE) {
  
  if (verbose) message("Visualizing ComplexMap with a color legend...")
  
  # Build base graph
  p <- .build_base_graph(layoutDf, edgesDf, edgeColor, verbose = FALSE)
  
  # Add node layer with appropriate color mapping
  if (is.null(color.by)) {
    p <- p + ggraph::geom_node_point(
      ggplot2::aes(fill = primaryFunctionalDomain, size = sizeMapping),
      shape = 21, color = edgeColor, stroke = 0.2
    )
  } else {
    p <- p + ggraph::geom_node_point(
      ggplot2::aes(fill = .data[[color.by]], size = sizeMapping),
      shape = 21, color = edgeColor, stroke = 0.2
    )
  }
  
  # Apply color scale
  p <- .apply_color_mapping(p, layoutDf, color.by, color.palette, 
                            color.legend.title, geom_type = "fill",
                            unenrichedColor = unenrichedColor)
  
  # Apply theme (with slight modification for fill legend)
  p <- p +
    ggplot2::scale_size_continuous(
      range = nodeSizeRange, name = size.legend.title
    ) +
    ggraph::scale_edge_alpha_continuous(
      range = c(0.1, 1), name = "Similarity"
    ) +
    ggraph::theme_graph(base_family = fontFamily) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bgColor),
      panel.background = ggplot2::element_rect(fill = bgColor),
      legend.background = ggplot2::element_rect(fill = bgColor),
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(color = edgeColor),
      legend.title = ggplot2::element_text(color = edgeColor),
      plot.title = ggplot2::element_text(color = edgeColor, hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(color = edgeColor, hjust = 0.5)
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)
  
  return(p)
}


#' Visualize a Complex Map Interactively
#'
#' @description
#' Creates a dynamic, zoomable HTML widget of the complex network using the
#' `visNetwork` package.
#'
#' @details
#' This function supports both categorical and continuous coloring. When `color.by`
#' is specified, nodes are colored by the numeric values in that column, and
#' the tooltip is updated to include this information.
#'
#' @inheritParams visualizeMapDirectLabels
#' @param width The width of the HTML widget.
#' @param height The height of the HTML widget.
#' @param physicsEnabled A logical value. If `TRUE`, nodes will physically
#'   react to dragging. Defaults to `FALSE` for a stable layout.
#' @param rescale_coords Logical. If `TRUE`, rescale coordinates to visNetwork's
#'   coordinate system. Defaults to `TRUE`.
#' @param rescale_range Numeric vector of length 2. Range for coordinate rescaling.
#'   Defaults to `c(-800, 800)`.
#' @param legend_steps Integer. Number of steps to show in continuous color legend.
#'   Defaults to 5.
#'
#' @return A `visNetwork` HTML widget object.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @export
#' @examples
#' # --- Sample Data ---
#' nodes <- tibble::tibble(
#'   complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
#'   primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Metabolism"),
#'   sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#0000FF"),
#'   proteinCount = c(10, 8, 12), proteins = c("A,B", "B,C", "D,E"),
#'   score = c(105.1, 88.3, 95.7)
#' )
#' edges <- tibble::tibble(
#'   source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
#' )
#'
#' # --- Generate Plot (requires visNetwork) ---
#' if (requireNamespace("visNetwork", quietly = TRUE)) {
#'   # Default categorical coloring
#'   visualizeMapInteractive(nodes, edges)
#'
#'   # Continuous coloring
#'   visualizeMapInteractive(nodes, edges, color.by = "score")
#' }
#'
visualizeMapInteractive <- function(
    layoutDf, edgesDf,
    width = "100%", height = "90vh",
    title = "ComplexMap Functional Landscape",
    physicsEnabled = FALSE,
    color.by = NULL, 
    color.palette = "viridis",
    rescale_coords = TRUE,
    rescale_range = c(-800, 800),
    legend_steps = 5,
    verbose = TRUE) {
  
  if (verbose) message("Generating interactive visNetwork plot...")
  
  # Prepare base node data
  nodesVn <- layoutDf %>%
    dplyr::mutate(
      id = complexId,
      label = complexId,
      value = sizeMapping
    )
  
  # Rescale coordinates if needed
  if (rescale_coords) {
    nodesVn <- nodesVn %>%
      dplyr::mutate(
        x = scales::rescale(x, to = rescale_range),
        y = scales::rescale(-y, to = rescale_range)
      )
  }
  
  # Apply color mapping and create tooltips
  if (is.null(color.by)) {
    # Categorical mode
    nodesVn <- nodesVn %>%
      dplyr::mutate(
        color = colorHex,
        title = .create_tooltip(complexId, primaryFunctionalDomain, 
                                proteinCount, proteins)
      )
  } else {
    # Continuous mode
    if (!color.by %in% names(layoutDf)) {
      stop(sprintf("Column '%s' not found in the provided node data frame.", color.by))
    }
    
    quant_values <- layoutDf[[color.by]]
    color_mapper <- scales::col_numeric(
      palette = color.palette,
      domain = range(quant_values, na.rm = TRUE)
    )
    
    nodesVn <- nodesVn %>%
      dplyr::mutate(
        color = ifelse(is.na(.data[[color.by]]), "#888888", 
                       color_mapper(.data[[color.by]])),
        title = .create_tooltip(complexId, primaryFunctionalDomain, 
                                proteinCount, proteins, 
                                color.by = color.by, 
                                color.value = .data[[color.by]])
      )
  }
  
  # Prepare edge data
  edgesVn <- edgesDf %>%
    dplyr::mutate(
      from = source_complex_id, 
      to = target_complex_id, 
      value = weight
    )
  
  # Create network object
  vn <- visNetwork::visNetwork(
    nodes = nodesVn, edges = edgesVn, 
    width = width, height = height, main = title
  ) %>%
    visNetwork::visNodes(shape = "dot", borderWidth = 2, shadow = TRUE) %>%
    visNetwork::visEdges(
      smooth = FALSE, color = list(color = "#888888", highlight = "#00BFFF")
    ) %>%
    visNetwork::visPhysics(enabled = physicsEnabled) %>%
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE,
      selectedBy = "primaryFunctionalDomain"
    ) %>%
    visNetwork::visLayout(randomSeed = 123) %>%
    visNetwork::visInteraction(
      navigationButtons = TRUE, dragNodes = TRUE, dragView = TRUE,
      zoomView = TRUE, tooltipDelay = 200
    )
  
  # Add continuous legend if needed
  if (!is.null(color.by)) {
    domain <- range(layoutDf[[color.by]], na.rm = TRUE)
    vals <- seq(domain[1], domain[2], length.out = legend_steps)
    color_mapper <- scales::col_numeric(color.palette, domain = domain)
    
    vn <- vn %>%
      visNetwork::visLegend(
        useGroups = FALSE, 
        addNodes = data.frame(
          label = round(vals, 2),
          color = color_mapper(vals),
          shape = "dot",
          size = 10
        ),
        main = list(
          text = color.by, 
          style = "color:black;font-weight:bold;font-size:14px;text-align:center;"
        )
      )
  }
  
  return(vn)
}