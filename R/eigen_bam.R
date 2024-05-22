#' eigen_bam: Compute the Eigen system of two bam objects
#' @description Calculates the Eigen values and Eigen vectors of bam objects
#'
#' @param A A bam object of class setA.
#' @param M A bam object of class setM.
#' @param AM A logical value to specify whether to use the product AM or MA.
#' If true the AM will be returned else the product MA will be returned.
#' @param which_eigen An integer representing the which eigen value and eigen
#' vector will be computed.
#' @param rmap Logical. If TRUE the function will return a map of the eigen
#' vector of the product AM.
#' @return A list with four objects. 1) eigen_values (these are indicated in
#' which_eigen parameter of the function), 2) eigen_vectors (the corresponding
#' eigen vectors of each eigen value), 3) Standardized eigen vectors (0 to 1),
#' 4) A RasterLayer depicting the information of the first eigen vector of
#' the system.
#' @details The eigenvector associated with the dominant eigenvalue of an
#' adjacency matrix provides information about the number of forms in
#' which a cell can be visited from other cells. Details about the
#' eigen analysis in the context of the area of distribution can be
#' found in Soberon and Osorio-Olvera (2022).
#' @references
#' \insertRef{SoberonOsorio}{bamm}.
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' \donttest{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#' sparse_mod <- bamm::model2sparse(model = model,0.75)
#' plot(sparse_mod@niche_model)
#' adj_mod <- bamm::adj_mat(sparse_mod,ngbs = 1,eigen_sys = TRUE)
#' # Product AM
#' eig_bam_am <- bamm::eigen_bam(A=sparse_mod,M=adj_mod,AM=TRUE)
#' raster::plot(eig_bam_am$map)
#' # Product MA
#' eig_bam_ma <- bamm::eigen_bam(A=sparse_mod,M=adj_mod,AM=FALSE)
#' raster::plot(eig_bam_ma$map)
#' }
#' @export

eigen_bam <- function(A=NULL,M=NULL,AM=TRUE, which_eigen=1,rmap=TRUE){
  if(!is.null(A) && !methods::is(A, "setA"))
    stop("Object A should be of class setA")
  if(!is.null(M) && !methods::is(M, "setM"))
    stop("Object M should be of class setM")

  normalize_vec <- function(x){
    x <- abs(x)
    return((x - min(x)) / (max(x) - min(x)))
  }

  if(AM){
    prod_mat <- A@sparse_model %*% M@adj_matrix
  } else{
    prod_mat <- M@adj_matrix %*% A@sparse_model
  }


  r1 <- RSpectra::eigs(prod_mat,which_eigen)

  eigen_vectors_norm <- apply(r1$vectors,
                              MARGIN = 2,
                              normalize_vec)

  r2 <- list(eigen_values = r1$values,
             eigen_vectors=r1$vectors,
             eigen_vectors_norm=eigen_vectors_norm)


  if(rmap){
    m1 <- A@niche_model
    m1[A@cellIDs] <- r2$eigen_vectors_norm[,which_eigen]
    r2$map <- m1
  }
  return(r2)
}
