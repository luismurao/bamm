#' permute_pam: Function to permute a Presence-Absence-Matrix.
#' @param m Presence-Absence-Matrix (PAM) or a binary matrix with columns
#'          representing species and rows sites.
#' @param niter Number of iterations to permute the PAM.
#' @param as_sparse If TRUE the PAM will be returned as a sparse matrix
#' @return Returns a permuted matrix of the same dimensions of m
#' (same number of rows and columns). Note that the sum of each row and column
#' of this permuted matrix is equal to that of m.
#' species.
#' @details This function is an implementation of the curve ball algorithm
#'           following Strona et al. (2014).
#' @references
#' \insertRef{Strona2014}{bamm}.
#'
#' @export
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' set.seed(111)
#' pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
#' ppam <- bamm::permute_pam(m = pam,niter = NULL,as_sparse = FALSE)
#' # Check if matrices are different
#' all(pam == ppam)
#' # Check if row totals are the same
#' all(Matrix::rowSums(pam) == Matrix::rowSums(ppam))
#' # Check if column total are the same
#' all(Matrix::colSums(pam) == Matrix::colSums(ppam))


permute_pam <- function(m,niter=NULL,as_sparse=FALSE){
  if(is.data.frame(m))
    m <- as.matrix(m)
  if(!is.matrix(m))
    stop("m should be a matrix or a data.frame")
  if(is.null(niter))
    niter <- nrow(m)*5
  if(!is.numeric(niter))
    stop("niter shuld be an integer")
  ppam <- permute_matrix(m,niter)
  if(as_sparse){
    ppam <- Matrix::Matrix(ppam, sparse = TRUE)
  }

  return(ppam)
}
