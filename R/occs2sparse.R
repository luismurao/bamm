#' occs2sparse: Converts occurrence data into a sparse matrix object
#' @param modelsparse A setA object returned by the function
#' \code{\link[bamm]{model2sparse}}
#' @param occs A matrix or a data.frame containing two columns.
#' The first one is the longitude and the second is the latitude.
#' @return A sparse vector of zeros (presences) and ones (absences).
#' @details  Rows of this column vector represent non NA pixels of the
#' niche model.
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @examples
#'
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#'
#' sparse_mod <- bamm::model2sparse(model,threshold=0.05)
#'
#' occs_lep_cal <- data.frame(longitude = c(-115.10417,
#'                                          -104.90417),
#'                            latitude = c(29.61846,
#'                                         29.81846))
#'
#' occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
#'                                 occs = occs_lep_cal)
#'
#' head(occs_sparse)

occs2sparse <- function(modelsparse,occs){
  if(!methods::is(modelsparse,"setA")){
    stop("modelsparse should be of class setA")
  }

  occsIDs <- raster::cellFromXY(modelsparse@niche_model,occs)
  occsIDsSparese <-  which(modelsparse@cellIDs %in% occsIDs)
  pres_abs <- numeric(nrow(modelsparse@sparse_model))
  pres_abs[occsIDsSparese] <- 1
  pres_abs <- Matrix::Matrix(data = pres_abs,sparse=TRUE)

  #mod_atts <- setA(niche_model = modelsparse@niche_model,
  #                 cellIDs = modelsparse@cellIDs,
  #                 sparse_model =  modelsparse@sparse_model,
  #                 occs_sparse = pres_abs,
  #                 n_occs = nrow(occs))
  return(pres_abs)
}
