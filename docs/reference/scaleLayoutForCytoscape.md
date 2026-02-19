# Scale Layout Coordinates for Cytoscape.js

Mean-centers and rescales raw layout coordinates to a fixed pixel range
suitable for Cytoscape.js.

## Usage

``` r
scaleLayoutForCytoscape(x, y, range = c(-500, 500))
```

## Arguments

- x:

  Numeric vector of x coordinates.

- y:

  Numeric vector of y coordinates.

- range:

  Numeric vector of length 2 giving target range. Default \`c(-500,
  500)\`.

## Value

A \`list\` with elements \`x\` and \`y\`, both rescaled.
