#' shape2Grid: Function to create a grid given a spatial polygon
#'
#' @description shapeToGrid creates a raster grid given a spatial polygon and a grid resolution.
#' @param shpolygon A SpatialPolygon, SpatialPolygonDataFrame representing the desired shape of the grid.
#' @param resolution Numeric. Spatial resoltion of the grid.
#' @param ones Logical. Fill with ones the values of the raster. If not the values will be written as cellID values.
#' @return Returns a raster object whith the shape of 'shpolygon' of a given resoltion.
#' @export
#' @examples
#' \dontrun{
#' data("wrld_simpl", package = "maptools")
#' mx <- wrld_simpl[wrld_simpl$NAME=="Mexico",]
#' mx_grid <- bamm::shape2Grid(mx,0.5)
#' }


shape2Grid <- function(shpolygon,resolution,ones=TRUE){
  r1 <- raster::raster()
  raster::extent(r1) <- raster::extent(shpolygon)
  raster::res(r1) <- resolution
  if(ones)
    r1[] <- 1
  else
    r1[] <- 1:raster::ncell(r1)
  r1 <- raster::crop(r1,shpolygon)
  r1 <- raster::mask(r1,shpolygon)
  return(r1)
}
