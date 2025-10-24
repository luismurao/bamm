#' Convert distribution polygons to a presence-absence matrix (PAM)
#'
#' @description
#' This function takes spatial polygon objects (typically species distribution polygons)
#' and converts them into a presence-absence matrix (PAM) using a rasterization approach
#' at a specified resolution. The function is particularly useful for biodiversity
#' and biogeography analyses.
#'
#' @param poly Spatial polygon object (SpatialPolygonsDataFrame, sf, etc.) containing
#'   species or taxon distribution data.
#' @param taxon_attribute Character. The name of the column in `poly` that contains
#'   taxon identifiers (species, genera, etc.).
#' @param resolution Numeric. The resolution for the output raster in coordinate
#'   system units (cell size).
#' @param polymask Optional. A spatial polygon object used to mask and crop the
#'   resulting raster. If NULL, no masking is applied.
#'
#' @return A PAM (Presence-Absence Matrix) object of class `pam` from the `bamm`
#'   package containing:
#'   \itemize{
#'     \item \code{matrix}: Presence-absence matrix (1/0) where rows represent grid cells
#'           and columns represent taxa
#'     \item \code{richness}: Richness pattern across grid cells
#'     \item \code{sparse}: Logical indicating if the matrix is stored in sparse format
#'     \item \code{cell_coordinates}: Coordinates of each grid cell
#'   }
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Creating a base raster with the specified resolution and extent
#'   \item Splitting polygons by taxon attribute
#'   \item Rasterizing each taxon's distribution using exact extraction
#'   \item Stacking individual rasters into a multi-layer raster
#'   \item Applying optional masking with polymask
#'   \item Converting the raster stack to a PAM using bamm::models2pam
#' }
#'
#' @note
#' This function requires the following packages: raster, exactextractr, purrr, and bamm.
#' The input polygons are converted to SpatialPolygonsDataFrame if they aren't already.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(raster)
#'
#' # Example with sample data
#' uicn <- readRDS(system.file("extdata/uicn.rds",package = "bamm"))
#' sudam <- readRDS(system.file("extdata/suam.rds",package = "bamm"))
#' # Convert to PAM with 0.5 degree resolution
#' pam_result <- bamm::pol2pam(poly = uicn,
#'                             taxon_attribute = "binomial",
#'                             resolution = 0.5,
#'                             polymask = NULL)
#'
#' # With masking polygon
#' pam_masked <- pol2pam(poly = uicn,
#'                       taxon_attribute = "binomial",
#'                       resolution = 0.5,
#'                       polymask = sudam)
#' }
#'
#' @seealso
#' \code{\link[bamm]{models2pam}}, \code{\link[raster]{raster}},
#' \code{\link[exactextractr]{exact_extract}}
#'
#' @export
#' @importFrom raster raster extent res stack mask crop
#' @importFrom exactextractr exact_extract
#' @importFrom purrr map compact
#' @importFrom methods as
pol2pam <- function(poly, taxon_attribute, resolution, polymask = NULL){
  r1 <- raster::raster()
  raster::extent(r1) <- raster::extent(poly)
  raster::res(r1) <- rep(resolution, 2)
  r1[] <- 0
  #poly <- poly |> as("Spatial")
  ex <- exactextractr::exact_extract(r1,
                               poly,
                               progress = FALSE,
                               include_cell = TRUE)
  names(ex) <- poly[[taxon_attribute]]
  bd <- seq_along(ex) |> purrr::map_dfr(function(x){
    if(nrow(ex[[x]])> 0L){
      data.frame(taxon = names(ex[x]),ex[[x]])
    }
  })
  bdL <- bd |> split.data.frame(bd$taxon)
  pol2ras <- seq_along(bdL) |> purrr::map(function(x){
    r1[bdL[[x]]$cell] <- 1
    names(r1) <- bdL[[x]]$taxon[1]

    return(r1)
  })
  pol2ras <- raster::stack(pol2ras)
  if(!is.null(polymask)){
    pol2ras <- raster::mask(raster::crop(pol2ras, polymask), polymask)
  }
  pol2pam <- bamm::models2pam(pol2ras, return_coords = TRUE, sparse = FALSE)
  return(pol2pam)
}
