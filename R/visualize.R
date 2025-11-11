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
#' This function is designed for clarity when domains are distinct. For large,
#' dense domains (more than 2 complexes), a single label is placed at the
#' domain's centroid. For smaller domains, each complex is labeled
#' individually. It uses the pre-computed layout from `computeMapTopology`.
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
#' edges <- tibble::tibble(
#'   source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
#' )
#'
#' # --- Generate Plot ---
#' if (requireNamespace("ggrepel", quietly = TRUE)) {
#'   visualizeMapDirectLabels(nodes, edges)
#' }
#'
#' @export
#'
visualizeMapDirectLabels <- function(
    layoutDf, edgesDf,
    title="ComplexMap Functional Landscape",
    subtitle="Nodes are protein complexes, colored by function",
    bgColor="black", edgeColor="white",
    nodeSizeRange=c(2, 10), labelFillColor=ggplot2::alpha("white", 0.7),
    fontFamily="sans", verbose=TRUE) {
  
  if (verbose) message("Visualizing ComplexMap with direct labels...")
  
  if (!requireNamespace("ggrepel", quietly=TRUE)) {
    stop("Package 'ggrepel' is required for direct labeling.", call.=FALSE)
  }
  
  graphForPlotting <- tidygraph::tbl_graph(
    nodes=layoutDf, edges=edgesDf, directed=FALSE
  )
  
  labelData <- layoutDf %>%
    dplyr::filter(
      primaryFunctionalDomain != "Unenriched" &
        !is.na(primaryFunctionalDomain)
    ) %>%
    dplyr::add_count(primaryFunctionalDomain, name="domainSize")
  
  centroidLabels <- labelData %>%
    dplyr::filter(domainSize > 2) %>%
    dplyr::group_by(primaryFunctionalDomain) %>%
    dplyr::summarise(x=mean(x), y=mean(y), .groups="drop")
  
  directLabels <- labelData %>% dplyr::filter(domainSize <= 2)
  
  ggraph::ggraph(graphForPlotting, layout='manual', x=x, y=y) +
    ggraph::geom_edge_fan(
      ggplot2::aes(alpha=weight), color=edgeColor, width=0.25
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size=sizeMapping, color=colorHex)
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_continuous(
      range=nodeSizeRange, name="Complex Size (log2)",
      guide=ggplot2::guide_legend(override.aes=list(color="white"))
    ) +
    ggraph::scale_edge_alpha_continuous(
      range=c(0.1, 1), name="Similarity"
    ) +
    ggrepel::geom_label_repel(
      data=centroidLabels,
      ggplot2::aes(x=x, y=y,
                   label=stringr::str_wrap(primaryFunctionalDomain, 20)),
      fill=labelFillColor, fontface='bold', color="black",
      segment.color=edgeColor
    ) +
    ggrepel::geom_label_repel(
      data=directLabels,
      ggplot2::aes(x=x, y=y,
                   label=stringr::str_wrap(primaryFunctionalDomain, 20)),
      fill=labelFillColor, fontface='bold', color="black",
      segment.color=edgeColor, box.padding=0.5, point.padding=0.5
    ) +
    ggraph::theme_graph(base_family=fontFamily) +
    ggplot2::theme(
      plot.background=ggplot2::element_rect(fill=bgColor),
      panel.background=ggplot2::element_rect(fill=bgColor),
      legend.background=ggplot2::element_rect(fill=bgColor),
      legend.key=ggplot2::element_rect(fill=bgColor),
      legend.text=ggplot2::element_text(color=edgeColor),
      legend.title=ggplot2::element_text(color=edgeColor),
      plot.title=ggplot2::element_text(color=edgeColor, hjust=0.5, size=16),
      plot.subtitle=ggplot2::element_text(color=edgeColor, hjust=0.5)
    ) +
    ggplot2::labs(title=title, subtitle=subtitle)
}

#' Visualize a Complex Map with a Color Legend
#'
#' @description
#' Creates a static visualization of the complex network using `ggraph`, where
#' functional domains are represented by a discrete color scale in a legend.
#'
#' @details
#' This plot is useful for overviews where direct labels would be too
#' cluttered. It maps the `primaryFunctionalDomain` to node color and displays
#' a legend. It uses the pre-computed layout from `computeMapTopology`.
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
#' edges <- tibble::tibble(
#'   source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
#' )
#'
#' # --- Generate Plot ---
#' visualizeMapWithLegend(nodes, edges)
#'
#' @export
#'
visualizeMapWithLegend <- function(
    layoutDf, edgesDf,
    title="ComplexMap Functional Landscape",
    subtitle="Nodes are protein complexes, colored by function",
    bgColor="black", edgeColor="white",
    nodeSizeRange=c(2, 10), unenrichedColor="#CCCCCC",
    fontFamily="sans", verbose=TRUE) {
  
  if (verbose) message("Visualizing ComplexMap with a color legend...")
  
  graphForPlotting <- tidygraph::tbl_graph(
    nodes=layoutDf, edges=edgesDf, directed=FALSE
  )
  
  domainColorMap <- layoutDf %>%
    dplyr::filter(
      primaryFunctionalDomain != "Unenriched" &
        !is.na(primaryFunctionalDomain)
    ) %>%
    dplyr::group_by(primaryFunctionalDomain) %>%
    dplyr::summarise(
      color=names(which.max(table(colorHex))), .groups="drop"
    )
  
  legendColors <- tibble::deframe(domainColorMap)
  legendColors["Unenriched"] <- unenrichedColor
  
  ggraph::ggraph(graphForPlotting, layout='manual', x=x, y=y) +
    ggraph::geom_edge_fan(
      ggplot2::aes(alpha=weight), color=edgeColor, width=0.25
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill=primaryFunctionalDomain, size=sizeMapping),
      shape=21, color=edgeColor, stroke=0.2
    ) +
    ggplot2::scale_fill_manual(
      name="Functional Domain", values=legendColors
    ) +
    ggplot2::scale_size_continuous(
      range=nodeSizeRange, name="Complex Size (log2)"
    ) +
    ggraph::scale_edge_alpha_continuous(
      range=c(0.1, 1), name="Similarity"
    ) +
    ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(size=5))) +
    ggraph::theme_graph(base_family=fontFamily) +
    ggplot2::theme(
      plot.background=ggplot2::element_rect(fill=bgColor),
      panel.background=ggplot2::element_rect(fill=bgColor),
      legend.background=ggplot2::element_rect(fill=bgColor),
      legend.key=ggplot2::element_blank(),
      legend.text=ggplot2::element_text(color=edgeColor),
      legend.title=ggplot2::element_text(color=edgeColor),
      plot.title=ggplot2::element_text(color=edgeColor, hjust=0.5, size=16),
      plot.subtitle=ggplot2::element_text(color=edgeColor, hjust=0.5)
    ) +
    ggplot2::labs(title=title, subtitle=subtitle)
}

#' Visualize a Complex Map Interactively
#'
#' @description
#' Creates a dynamic, zoomable, and explorable HTML widget of the complex
#' network using the `visNetwork` package.
#'
#' @details
#' This function is ideal for data exploration. Nodes can be clicked and
#' dragged, and hovering over a node reveals a detailed tooltip with its
#' attributes. The coordinates from the static layout are rescaled and used
#' to provide an initial, non-random arrangement.
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
#'
visualizeMapInteractive <- function(
    layoutDf, edgesDf,
    width="100%", height="90vh",
    title="ComplexMap Functional Landscape",
    physicsEnabled=FALSE, verbose=TRUE) {
  
  if (verbose) message("Generating interactive visNetwork plot...")
  
  nodesVn <- layoutDf %>%
    dplyr::mutate(
      id = complexId,
      label = complexId,
      value = sizeMapping,
      color = colorHex,
      x = scales::rescale(x, to=c(-800, 800)),
      y = scales::rescale(-y, to=c(-800, 800)),
      title = paste0(
        "<div style='font-family:sans-serif; text-align:left;'>",
        "<b>Complex:</b> ", complexId, "<br>",
        "<b>Function:</b> ", primaryFunctionalDomain, "<br>",
        "<b>Protein Count:</b> ", proteinCount, "<hr>",
        "<b>Members:</b><br>", gsub(",", ", ", proteins), "</div>"
      )
    )
  
  edgesVn <- edgesDf %>%
    dplyr::mutate(
      from=source_complex_id, to=target_complex_id, value=weight
    )
  
  visNetwork::visNetwork(
    nodes=nodesVn, edges=edgesVn, width=width, height=height, main=title
  ) %>%
    visNetwork::visNodes(shape="dot", borderWidth=2, shadow=TRUE) %>%
    visNetwork::visEdges(
      smooth=FALSE, color=list(color="#888888", highlight="#00BFFF")
    ) %>%
    visNetwork::visPhysics(enabled=physicsEnabled) %>%
    visNetwork::visOptions(
      highlightNearest=list(enabled=TRUE, degree=1, hover=TRUE),
      nodesIdSelection=TRUE,
      selectedBy="primaryFunctionalDomain"
    ) %>%
    visNetwork::visLayout(randomSeed=123) %>%
    visNetwork::visInteraction(
      navigationButtons=TRUE, dragNodes=TRUE, dragView=TRUE,
      zoomView=TRUE, tooltipDelay=200
    )
}