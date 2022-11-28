#' bam_islands: Estimates the number of connected patches of a distribution model

#' @description Estimates the number of connected patches of a distribution model
#' @param model A niche model projected in an area of interest (could be the M).
#' @param ngbs Number of pixels the species can disparse.
#' @export
#' @examples

#'
#' \dontrun{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#' model <- model > 0.7
#' islands <- bamm::bam_islands(model,ngbs=1)
#' }
#' @return A list with the product of AM (A abiotic niche, M accesible area);
#' the degree matrix of the AM matrix (see)
#'
bam_islands <- function(model,ngbs=1){
  sparse_mod <- bamm::model2sparse(model)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs=ngbs)
  ama <- sparse_mod@sparse_model %*% adj_mod@adj_matrix %*%
    sparse_mod@sparse_model

  nr <- nrow(ama)

  grado <- Matrix::rowSums(ama)
  m_grado <- sparse_mod@sparse_model
  Matrix::diag(m_grado)  <- grado
  Lp <- m_grado - ama
  message(paste("Computing Eigen system for Laplacian matrix" ,
                "<this process could take time>... "))
  eigsLp <- RSpectra::eigs(Lp,ncol(ama))
  no_nicho <- length(which(Matrix::diag(sparse_mod@sparse_model) == 0))

  eigVals <- round(abs(eigsLp$values),5)
  evec_ceros <-  length(which(eigVals==0))
  islas <- evec_ceros - no_nicho

  results <- list(ama=ama,
                  degree=m_grado,
                  laplacian_matrix=Lp,
                  eigs_sys_lp = eigsLp,
                  islands=islas)
  return(results)

}
