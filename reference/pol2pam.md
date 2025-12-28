# Convert distribution polygons to a presence-absence matrix (PAM)

This function takes spatial polygon objects (typically species
distribution polygons) and converts them into a presence-absence matrix
(PAM) using a rasterization approach at a specified resolution. The
function is particularly useful for biodiversity and biogeography
analyses.

## Usage

``` r
pol2pam(poly, taxon_attribute, resolution, polymask = NULL)
```

## Arguments

- poly:

  Spatial polygon object (SpatialPolygonsDataFrame, sf, etc.) containing
  species or taxon distribution data.

- taxon_attribute:

  Character. The name of the column in \`poly\` that contains taxon
  identifiers (species, genera, etc.).

- resolution:

  Numeric. The resolution for the output raster in coordinate system
  units (cell size).

- polymask:

  Optional. A spatial polygon object used to mask and crop the resulting
  raster. If NULL, no masking is applied.

## Value

A PAM (Presence-Absence Matrix) object of class \`pam\` from the
\`bamm\` package containing:

- `matrix`: Presence-absence matrix (1/0) where rows represent grid
  cells and columns represent taxa

- `richness`: Richness pattern across grid cells

- `sparse`: Logical indicating if the matrix is stored in sparse format

- `cell_coordinates`: Coordinates of each grid cell

## Details

The function works by:

1.  Creating a base raster with the specified resolution and extent

2.  Splitting polygons by taxon attribute

3.  Rasterizing each taxon's distribution using exact extraction

4.  Stacking individual rasters into a multi-layer raster

5.  Applying optional masking with polymask

6.  Converting the raster stack to a PAM using bamm::models2pam

## Note

This function requires the following packages: raster, exactextractr,
purrr, and bamm. The input polygons are converted to
SpatialPolygonsDataFrame if they aren't already.

## See also

[`models2pam`](models2pam.md),
[`raster`](https://rdrr.io/pkg/raster/man/raster.html),
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load required libraries
library(raster)

# Example with sample data
uicn <- readRDS(system.file("extdata/uicn.rds",package = "bamm"))
sudam <- readRDS(system.file("extdata/suam.rds",package = "bamm"))
# Convert to PAM with 0.5 degree resolution
pam_result <- bamm::pol2pam(poly = uicn,
                            taxon_attribute = "binomial",
                            resolution = 0.5,
                            polymask = NULL)

# With masking polygon
pam_masked <- pol2pam(poly = uicn,
                      taxon_attribute = "binomial",
                      resolution = 0.5,
                      polymask = sudam)
} # }
```
