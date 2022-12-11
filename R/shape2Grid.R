#' shape2Grid: Function to create a grid given a spatial polygon
#'
#' @description shape2Grid creates a raster grid given a spatial polygon and
#' a grid resolution.
#' @param shpolygon A SpatialPolygon, SpatialPolygonDataFrame representing
#' the desired shape of the grid.
#' @param resolution Numeric. Spatial resolution of the grid.
#' @param ones Logical. Fill with ones the values of the raster. If not the
#' values will be written as cellID values.
#' @return Returns a raster object with the shape of 'shpolygon' of a given
#'  resolution.
#' @export
#' @importFrom raster raster
#' @importFrom raster extent
#' @importFrom raster res
#' @importFrom raster crop
#' @importFrom raster mask
#' @importFrom sp Polygon
#' @examples
#' \dontrun{
#' x_coord <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
#' y_coord <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125,
#'              59.736328125)
#' xy <- cbind(x_coord, y_coord)
#' p <- sp::Polygon(xy)
#' ps <- sp::Polygons(list(p),1)
#' sps <- sp::SpatialPolygons(list(ps))
#' r1 <- bamm::shape2Grid(sps,resolution = 0.5,ones = T)
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
