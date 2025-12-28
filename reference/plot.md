# Plot method for objects of class diversity_range bamm.

Plot method for objects of class diversity_range bamm.

## Usage

``` r
# S4 method for class 'diversity_range,ANY'
plot(
  x,
  xlab = NULL,
  plot_type = "diversity_range",
  legend = TRUE,
  legend_position = "bottomright",
  ylab = NULL,
  col = NULL,
  pch = NULL,
  pch_legend = 19,
  radius = 0.5,
  ...
)
```

## Arguments

- x:

  An object of class diversity_range

- xlab:

  x label

- plot_type:

  Plot type: possible options: "diversity_range" (range-diversity plot),
  "diversity_range_map" (a raster map with diversity_range categories),
  "alpha" (a raster map with alpha diversity values), "dispersion_field"
  (a raster with dispersion field)

- legend:

  Logical. If TRUE the legend of the categorical diversity range values
  will appear.

- legend_position:

  Legend position.

- ylab:

  y label

- col:

  Plot colors.

- pch:

  Patch type.

- pch_legend:

  Patch type for legends.

- radius:

  Size of the patch for the interactive map.

- ...:

  Graphical parameters. Any argument that can be passed to 1)
  base::plot, such as axes=FALSE, main='title', ylab='latitude' 2)
  leaflet::leaflet or 3) leaflet::addCircleMarkers.

## Value

Plot of the results of the diversity_range analysis

## Details

To show interactive diversity_range plots install the 'plotly' R
package.

## Author

Luis Osorio-Olvera & Jorge Soberón
