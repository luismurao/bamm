#' permute_pam: Function to permute a Presence-Absence-Matrix.
#' @param m Presence-Absence-Matrix (PAM) or a binary matrix with columns
#'          representing species and rows sites.
#' @param niter Number of iterations to permute the PAM.
#' @param as_sparse If TRUE the PAM will be returned as a sparse matrix
#' @param randal Randomization algorithm applied to the PAM.
#' Possible choices "curveball", "fastball", and "indep_swap".
#' @return Returns a permuted matrix of the same dimensions of m
#' (same number of rows and columns). Note that the sum of each row and column
#' of this permuted matrix is equal to that of m.
#' species.
#' @details This function can use the "curveball" (Strona et al., 2014), the
#'          fastball (Godard and Neal, 2022), and the independent swap
#'          algorithms. The implementation of the "fastball" in C++ is provided
#'          in \url{https://github.com/zpneal/fastball/blob/main/fastball.cpp}.
#'          Please when using the "fastball" algorithm for publications cite
#'          Godard and Neal (2022). When using the "curveball" cite
#'          Strona et al. (2014). When using independent swap ("indep_swap")
#'          cite Kembel et al. (2010)
#'
#' @references
#' \insertRef{Strona2014}{bamm}
#'
#' \insertRef{Gordard2022}{bamm}
#'
#' \insertRef{Kembel2010}{bamm}
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


permute_pam <- function(m,niter=NULL,as_sparse=FALSE,randal="indep_swap"){
  ral <- match.arg(arg = randal,
                   choices = c("indep_swap","curveball","fastball"))
  if(is.data.frame(m))
    m <- as.matrix(m)
  if(!methods::is(m,"matrix") & !is.numeric(m[1,1])){
    stop("pam object should be a binary matrix")
  }
  if(is.null(niter))
    niter <- nrow(m)*5
  if(!is.numeric(niter))
    stop("niter shuld be an integer")
  if(randal == "curveball"){
    ppam <- permute_matrix(m,niter)
  } else if(randal == "fastball"){
    ppam <- permute_matrix_fb(m,niter)
  } else if(randal == "indep_swap"){
    ppam <- permute_matrix_indswap(m,niter)
  }#else{
    #warning("Algorithm not available:\n running the fastball algorithm...",
    #        call. = TRUE)
    #ppam <- permute_matrix_fb(m,niter)
  #}
  if(as_sparse){
    ppam <- Matrix::Matrix(ppam, sparse = TRUE)
  }

  return(ppam)
}
