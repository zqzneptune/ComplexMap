# Visualize a Complex Map Interactively

Creates a dynamic, zoomable HTML widget of the complex network using the
\`visNetwork\` package.

## Usage

``` r
visualizeMapInteractive(
  layoutDf,
  edgesDf,
  width = "100%",
  height = "90vh",
  title = "ComplexMap Functional Landscape",
  physicsEnabled = FALSE,
  color.by = NULL,
  color.palette = "viridis",
  rescale_coords = TRUE,
  rescale_range = c(-800, 800),
  legend_steps = 5,
  verbose = TRUE
)
```

## Arguments

- layoutDf:

  A data frame containing node attributes and layout coordinates,
  typically from \`computeMapTopology\`.

- edgesDf:

  A data frame containing the network edges.

- width:

  The width of the HTML widget.

- height:

  The height of the HTML widget.

- title:

  The main title for the plot.

- physicsEnabled:

  A logical value. If \`TRUE\`, nodes will physically react to dragging.
  Defaults to \`FALSE\` for a stable layout.

- color.by:

  A character string specifying the name of a numeric column in
  \`layoutDf\` to use for continuous node coloring. If \`NULL\`
  (default), the categorical \`colorHex\` column is used.

- color.palette:

  A character string or vector of colors for the continuous gradient
  (e.g., "viridis", "plasma", or \`c("blue", "white", "red")\`). Only
  used when \`color.by\` is specified.

- rescale_coords:

  Logical. If \`TRUE\`, rescale coordinates to visNetwork's coordinate
  system. Defaults to \`TRUE\`.

- rescale_range:

  Numeric vector of length 2. Range for coordinate rescaling. Defaults
  to \`c(-800, 800)\`.

- legend_steps:

  Integer. Number of steps to show in continuous color legend. Defaults
  to 5.

- verbose:

  A logical value indicating whether to print progress messages.

## Value

A \`visNetwork\` HTML widget object.

## Details

This function supports both categorical and continuous coloring. When
\`color.by\` is specified, nodes are colored by the numeric values in
that column, and the tooltip is updated to include this information.

## Author

Qingzhou Zhang \<zqzneptune@hotmail.com\>

## Examples

``` r
# --- Sample Data ---
nodes <- tibble::tibble(
  complexId = c("C1", "C2", "C3"), x = c(1, 2, 1.5), y = c(1, 1, 2),
  primaryFunctionalDomain = c("DNA Repair", "DNA Repair", "Metabolism"),
  sizeMapping = c(3, 4, 3.5), colorHex = c("#FF0000", "#FF0000", "#0000FF"),
  proteinCount = c(10, 8, 12), proteins = c("A,B", "B,C", "D,E"),
  score = c(105.1, 88.3, 95.7)
)
edges <- tibble::tibble(
  source_complex_id = "C1", target_complex_id = "C2", weight = 0.8
)

# --- Generate Plot (requires visNetwork) ---
if (requireNamespace("visNetwork", quietly = TRUE)) {
  # Default categorical coloring
  visualizeMapInteractive(nodes, edges)

  # Continuous coloring
  visualizeMapInteractive(nodes, edges, color.by = "score")
}
#> Generating interactive visNetwork plot...
#> Generating interactive visNetwork plot...

{"x":{"nodes":{"complexId":["C1","C2","C3"],"x":[-800,800,0],"y":[800,800,-800],"primaryFunctionalDomain":["DNA Repair","DNA Repair","Metabolism"],"sizeMapping":[3,4,3.5],"colorHex":["#FF0000","#FF0000","#0000FF"],"proteinCount":[10,8,12],"proteins":["A,B","B,C","D,E"],"score":[105.1,88.3,95.7],"id":["C1","C2","C3"],"label":["C1","C2","C3"],"value":[3,4,3.5],"color":["#FDE725","#440154","#26828E"],"title":["<div style='font-family:sans-serif; text-align:left;'><b>Complex:<\/b> C1<br><b>Function:<\/b> DNA Repair<br><b>score:<\/b> 105.1<br><b>Protein Count:<\/b> 10<hr><b>Members:<\/b><br>A, B<\/div>","<div style='font-family:sans-serif; text-align:left;'><b>Complex:<\/b> C2<br><b>Function:<\/b> DNA Repair<br><b>score:<\/b> 88.3<br><b>Protein Count:<\/b> 8<hr><b>Members:<\/b><br>B, C<\/div>","<div style='font-family:sans-serif; text-align:left;'><b>Complex:<\/b> C3<br><b>Function:<\/b> Metabolism<br><b>score:<\/b> 95.7<br><b>Protein Count:<\/b> 12<hr><b>Members:<\/b><br>D, E<\/div>"]},"edges":{"source_complex_id":["C1"],"target_complex_id":["C2"],"weight":[0.8],"from":["C1"],"to":["C2"],"value":[0.8]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot","borderWidth":2,"shadow":true},"manipulation":{"enabled":false},"edges":{"color":{"color":"#888888","highlight":"#00BFFF"},"smooth":false},"physics":{"enabled":false},"interaction":{"hover":true,"dragNodes":true,"dragView":true,"navigationButtons":true,"tooltipDelay":200,"zoomView":true,"zoomSpeed":1},"layout":{"randomSeed":123}},"groups":null,"width":"100%","height":"90vh","idselection":{"enabled":true,"style":"width: 150px; height: 26px","useLabels":true,"main":"Select by id"},"byselection":{"enabled":true,"style":"width: 150px; height: 26px","multiple":false,"hideColor":"rgba(200,200,200,0.5)","highlight":false,"variable":"primaryFunctionalDomain","main":"Select by primaryFunctionalDomain","values":["DNA Repair","Metabolism"]},"main":{"text":"ComplexMap Functional Landscape","style":"font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:20px;text-align:center;"},"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","highlight":{"enabled":true,"hoverNearest":true,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);","legend":{"width":0.2,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"nodes":{"label":[88.3,92.5,96.7,100.9,105.1],"color":["#440154","#3B528B","#21908D","#5DC863","#FDE725"],"shape":["dot","dot","dot","dot","dot"],"size":[10,10,10,10,10]},"nodesToDataframe":true,"main":{"text":"score","style":"color:black;font-weight:bold;font-size:14px;text-align:center;"}}},"evals":[],"jsHooks":[]}
```
