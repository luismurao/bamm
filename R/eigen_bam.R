#' eigen_bam: Compute the Eigen system of two bam objects
#' @description Calucultates the Eigen values and Eigen vectors of bam objects
#'
#' @param A A bam object of class setA.
#' @param M A bam object of class setM.
#' @param which_eigen An integer representing the which eigen value and eigen vector will be computed.
#' @param rmap Logical. If TRUE the function will return a map of the eigen vector of the product AM.
#' @return A table
#' @examples
#' \dontrun{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#' sparse_mod <- bamm::model2sparse(model = model)
#' adj_mod <- bamm::adj_mat(sparse_mod,nbgs=1,eigen_sys = T)
#' eig_bam <- bamm::eigen_bam(A=sparse_mod,M=adj_mod)
#' raster::plot(eig_bam$map)
#' }
#' @export

eigen_bam <- function(A=NULL,M=NULL,which_eigen=1,rmap=TRUE){
  if(!is.null(A) && !methods::is(A, "setA"))
    stop("Object A should be of class setA")
  if(!is.null(M) && !methods::is(M, "setM"))
    stop("Object M should be of class setM")

  normalize_vec <- function(x){
    x <- abs(x)
    return((x - min(x)) / max(x) - min(x))
  }
  if(exists("A") && exists("M")){
    AM <- A@sparse_model %*% M@adj_matrix

    r1 <- RSpectra::eigs(AM,which_eigen)

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
  return()
}
