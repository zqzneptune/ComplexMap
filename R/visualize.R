utils::globalVariables(c("primaryFunctionalDomain", "domainSize", "x", "y", 
                         "weight", "sizeMapping", "colorHex", "proteins", 
                         "source_complex_id", "target_complex_id"))
#' Visualize a Complex Map with Direct Node Labels
#'
#' @description
#' Creates a static visualization of the complex network using `ggraph`, where
#' functional domain labels are placed directly on the plot near the nodes.
#'
#' @details
#' This function is designed for clarity when domains are distinct. For large, dense
#' domains (more than 2 complexes), a single label is placed at the domain's
#' centroid. For smaller domains, each complex is labeled individually. It uses
#' the pre-computed layout from `computeMapTopology`.
#'
#' @param layoutDf A data frame containing node attributes and layout coordinates,
#'   typically from `computeMapTopology`.
#' @param edgesDf A data frame containing the network edges.
#' @param title The main title for the plot.
#' @param subtitle The subtitle for the plot.
#' @param bgColor The background color of the plot.
#' @param edgeColor The color of the network edges.
#' @param nodeSizeRange A numeric vector of length 2 specifying the min and max
#'   node size.
#' @param labelFillColor The background fill color for the labels.
#' @param fontFamily The base font family for all plot text.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A `ggplot` object representing the network visualization.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @examples
#' # --- Sample Data ---
#' nodes <- tibble::tibble(
#'   complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
#'   primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Metabolism"),
#'   sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#0000FF")
#' )
#' edges <- tibble::tibble(source_complex_id = "C1", target_complex_id = "C2", weight = 0.8)
#'
#' # --- Generate Plot ---
#' if (requireNamespace("ggrepel", quietly = TRUE)) {
#'   visualizeMapDirectLabels(nodes, edges)
#' }
#'
#' @export
#' @importFrom tidygraph tbl_graph
#' @importFrom dplyr filter add_count group_by summarise
#' @importFrom ggraph ggraph create_layout geom_edge_fan geom_node_point theme_graph scale_edge_alpha_continuous
#' @importFrom ggplot2 aes scale_color_identity scale_size_continuous guide_legend
#' @importFrom ggplot2 alpha theme element_rect
#' @importFrom ggplot2 element_text labs
#' @importFrom stringr str_wrap
#'
visualizeMapDirectLabels <- function(
    layoutDf, edgesDf,
    title = "ComplexMap Functional Landscape",
    subtitle = "Nodes are protein complexes, colored by function",
    bgColor = "black", edgeColor = "white",
    nodeSizeRange = c(2, 10), labelFillColor = alpha("white", 0.7),
    fontFamily = "sans", verbose = TRUE) {
  
  if (verbose) message("Visualizing ComplexMap with direct labels...")
  
  # The ggrepel package is essential for this plot
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for direct labeling.", call. = FALSE)
  }
  
  graphForPlotting <- tidygraph::tbl_graph(
    nodes = layoutDf, edges = edgesDf, directed = FALSE
  )
  
  # Prepare data for labeling: separate large domain centroids from small ones
  labelData <- layoutDf %>%
    filter(
      primaryFunctionalDomain != "Unenriched" &
        !is.na(primaryFunctionalDomain)
    ) %>%
    add_count(primaryFunctionalDomain, name = "domainSize")
  
  centroidLabels <- labelData %>%
    filter(domainSize > 2) %>%
    group_by(primaryFunctionalDomain) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop")
  
  directLabels <- labelData %>% filter(domainSize <= 2)
  
  ggraph(graphForPlotting, layout = 'manual', x = x, y = y) +
    geom_edge_fan(aes(alpha = weight), color = edgeColor, width = 0.25) +
    geom_node_point(aes(size = sizeMapping, color = colorHex)) +
    scale_color_identity() +
    scale_size_continuous(range = nodeSizeRange, name = "Complex Size (log2)",
                          guide = guide_legend(override.aes = list(color = "white"))) +
    scale_edge_alpha_continuous(range = c(0.1, 1), name = "Similarity") +
    ggrepel::geom_label_repel(data = centroidLabels,
                              aes(x = x, y = y, label = str_wrap(primaryFunctionalDomain, 20)),
                              fill = labelFillColor, fontface = 'bold', color = "black",
                              segment.color = edgeColor) +
    ggrepel::geom_label_repel(data = directLabels,
                              aes(x = x, y = y, label = str_wrap(primaryFunctionalDomain, 20)),
                              fill = labelFillColor, fontface = 'bold', color = "black",
                              segment.color = edgeColor, box.padding = 0.5, point.padding = 0.5) +
    theme_graph(base_family = fontFamily) +
    theme(
      plot.background = element_rect(fill = bgColor),
      panel.background = element_rect(fill = bgColor),
      legend.background = element_rect(fill = bgColor),
      legend.key = element_rect(fill = bgColor),
      legend.text = element_text(color = edgeColor),
      legend.title = element_text(color = edgeColor),
      plot.title = element_text(color = edgeColor, hjust = 0.5, size = 16),
      plot.subtitle = element_text(color = edgeColor, hjust = 0.5)
    ) +
    labs(title = title, subtitle = subtitle)
}


#' Visualize a Complex Map with a Color Legend
#'
#' @description
#' Creates a static visualization of the complex network using `ggraph`, where
#' functional domains are represented by a discrete color scale in a legend.
#'
#' @details
#' This plot is useful for overviews where direct labels would be too cluttered.
#' It maps the `primaryFunctionalDomain` to node color and displays a legend.
#' It uses the pre-computed layout from `computeMapTopology`.
#'
#' @inheritParams visualizeMapDirectLabels
#' @param unenrichedColor The color for nodes in the "Unenriched" category.
#'
#' @return A `ggplot` object representing the network visualization.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @examples
#' # --- Sample Data ---
#' nodes <- tibble::tibble(
#'   complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
#'   primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Unenriched"),
#'   sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#CCCCCC")
#' )
#' edges <- tibble::tibble(source_complex_id = "C1", target_complex_id = "C2", weight = 0.8)
#'
#' # --- Generate Plot ---
#' visualizeMapWithLegend(nodes, edges)
#'
#' @export
#' @importFrom tibble deframe
#' @importFrom ggplot2 scale_color_manual guides
#'
visualizeMapWithLegend <- function(
    layoutDf, edgesDf,
    title = "ComplexMap Functional Landscape",
    subtitle = "Nodes are protein complexes, colored by function",
    bgColor = "black", edgeColor = "white",
    nodeSizeRange = c(2, 10), unenrichedColor = "#CCCCCC",
    fontFamily = "sans", verbose = TRUE) {
  
  if (verbose) message("Visualizing ComplexMap with a color legend...")
  
  graphForPlotting <- tidygraph::tbl_graph(
    nodes = layoutDf, edges = edgesDf, directed = FALSE
  )
  
  # Create a named vector for the color scale
  domainColorMap <- layoutDf %>%
    filter(
      primaryFunctionalDomain != "Unenriched" &
        !is.na(primaryFunctionalDomain)
    ) %>%
    group_by(primaryFunctionalDomain) %>%
    summarise(color = names(which.max(table(colorHex))), .groups = "drop")
  
  legendColors <- deframe(domainColorMap)
  legendColors["Unenriched"] <- unenrichedColor
  
  ggraph(graphForPlotting, layout = 'manual', x = x, y = y) +
    geom_edge_fan(aes(alpha = weight), color = edgeColor, width = 0.25) +
    
    # --- START OF FIX ---
    # 1. Use shape=21, which has a separate fill and color (outline)
    # 2. Map the domain to the `fill` aesthetic
    # 3. Set the outline `color` to a constant visible color
    geom_node_point(aes(fill = primaryFunctionalDomain, size = sizeMapping),
                    shape = 21, color = edgeColor, stroke = 0.2) +
    
    # 4. Use scale_fill_manual to apply the colors to the fill
    scale_fill_manual(name = "Functional Domain", values = legendColors) +
    # --- END OF FIX ---
    
    scale_size_continuous(range = nodeSizeRange, name = "Complex Size (log2)") +
    scale_edge_alpha_continuous(range = c(0.1, 1), name = "Similarity") +
    guides(fill = guide_legend(override.aes = list(size = 5))) + # Guide now refers to `fill`
    theme_graph(base_family = fontFamily) +
    theme(
      plot.background = element_rect(fill = bgColor),
      panel.background = element_rect(fill = bgColor),
      legend.background = element_rect(fill = bgColor),
      legend.key = element_blank(),
      legend.text = element_text(color = edgeColor),
      legend.title = element_text(color = edgeColor),
      plot.title = element_text(color = edgeColor, hjust = 0.5, size = 16),
      plot.subtitle = element_text(color = edgeColor, hjust = 0.5)
    ) +
    labs(title = title, subtitle = subtitle)
}

#' Visualize a Complex Map Interactively
#'
#' @description
#' Creates a dynamic, zoomable, and explorable HTML widget of the complex
#' network using the `visNetwork` package.
#'
#' @details
#' This function is ideal for data exploration. Nodes can be clicked and dragged,
#' and hovering over a node reveals a detailed tooltip with its attributes. The
#' coordinates from the static layout are rescaled and used to provide an
#' initial, non-random arrangement.
#'
#' @param layoutDf A data frame of node attributes from `computeMapTopology`.
#' @param edgesDf A data frame of network edges.
#' @param width The width of the HTML widget.
#' @param height The height of the HTML widget.
#' @param title The main title of the network visualization.
#' @param physicsEnabled A logical value. If `TRUE`, nodes will physically
#'   react to dragging. Defaults to `FALSE` for a stable layout.
#' @param verbose A logical value indicating whether to print progress messages.
#'
#' @return A `visNetwork` HTML widget object.
#'
#' @author Qingzhou Zhang <zqzneptune@hotmail.com>
#'
#' @examples
#' # --- Sample Data ---
#' nodes <- tibble::tibble(
#'   complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
#'   primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Metabolism"),
#'   sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#0000FF"),
#'   proteinCount = c(10, 8, 12), proteins = c("A,B", "B,C", "D,E")
#' )
#' edges <- tibble::tibble(
#'   source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
#' )
#'
#' # --- Generate Plot ---
#' if (requireNamespace("visNetwork", quietly = TRUE)) {
#'   visualizeMapInteractive(nodes, edges)
#' }
#'
#' @export
#' @importFrom visNetwork visNetwork visNodes visEdges visPhysics visOptions
#' @importFrom visNetwork visLayout visInteraction
#' @importFrom scales rescale
#'
visualizeMapInteractive <- function(
    layoutDf, edgesDf,
    width = "100%", height = "90vh",
    title = "ComplexMap Functional Landscape",
    physicsEnabled = FALSE, verbose = TRUE) {
  
  if (verbose) message("Generating interactive visNetwork plot...")
  
  # Prepare node data for visNetwork, creating HTML tooltips
  nodesVn <- layoutDf %>%
    mutate(
      id = complexId, label = complexId, value = sizeMapping,
      color = colorHex,
      x = rescale(x, to = c(-800, 800)),
      y = rescale(-y, to = c(-800, 800)),
      title = paste0(
        "<div style='font-family:sans-serif; text-align:left;'>",
        "<b>Complex:</b> ", complexId, "<br>",
        "<b>Function:</b> ", primaryFunctionalDomain, "<br>",
        "<b>Protein Count:</b> ", proteinCount, "<hr>",
        "<b>Members:</b><br>", gsub(",", ", ", proteins), "</div>"
      )
    )
  
  # Prepare edge data
  edgesVn <- edgesDf %>%
    mutate(from = source_complex_id, to = target_complex_id, value = weight)
  
  visNetwork(nodes = nodesVn, edges = edgesVn, width = width, height = height,
             main = title) %>%
    visNodes(shape = "dot", borderWidth = 2, shadow = TRUE) %>%
    visEdges(smooth = FALSE, color = list(color = "#888888", highlight = "#00BFFF")) %>%
    visPhysics(enabled = physicsEnabled) %>%
    visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE,
      selectedBy = "primaryFunctionalDomain"
    ) %>%
    visLayout(randomSeed = 123) %>%
    visInteraction(
      navigationButtons = TRUE, dragNodes = TRUE, dragView = TRUE,
      zoomView = TRUE, tooltipDelay = 200
    )
}