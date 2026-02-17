#' model2sparse: Converts a niche model into a diagonal sparse matrix
#' @param model A raster object representing the geographic projection
#' of a niche model.
#' @param threshold A threshold to convert a continuous model into a
#' binary model.
#' @importFrom Matrix sparseMatrix
#' @return An object of class \code{\link[bamm]{setA}}. The niche model
#' is stored as diagonal sparse matrix (slot sparse_model).
#' @details threshold parameter represents the suitability value used to
#' convert continuous model into a binary model.
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @examples
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#'
#' sparse_mod <- bamm::model2sparse(model, threshold=0.75)
#' print(sparse_mod)
#' raster::plot(sparse_mod@niche_model)
model2sparse <- function(model, threshold=NULL){

  is_continous <- function(model){
    ff <- raster::freq(model,digits=2, useNA="no")
    if(nrow(ff)>3)
      return(TRUE)
    else
      return(FALSE)
  }

  source_model <- model

  if(is.numeric(threshold)){
    model <- model >= threshold
    source_model <- source_model*model
  }

  is_cont <- is_continous(model = model)
  if(is_cont){
    base::stop("Please provide a suitability value to binarize model")
  }
  model_vals <- raster::getValues(model)*1
  in_calArea <- which(!is.na(model_vals))
  ncols <- nrows <- length(in_calArea)
  mod_sparse <- Matrix::sparseMatrix(i=seq_along(in_calArea),
                                     j=seq_along(in_calArea),
                                     x=model_vals[in_calArea],
                                     dims = c(nrows,ncols))*1
  mod_coords <- raster::coordinates(model)[in_calArea,]
  mod_atts <- setA(niche_model = model,
                   suit_threshold = ifelse(is.numeric(threshold),
                                           threshold,as.numeric(NA)),
                   cellIDs = in_calArea,
                   suit_values = if(is_continous(source_model)){
                     source_model[in_calArea]} else as.numeric(NA),
                   sparse_model =  mod_sparse,
                   coordinates= mod_coords)
  return(mod_atts)

}
