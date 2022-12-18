#' null_dispersion_field_distribution: Null distribution of the
#' dispersion field
#' @description null_dispersion_field_distribution estimates a
#' random distribution of the dispersion field values.
#' @param pam A Presence-Absence-Matrix of matrix class or sparse matrix.
#' @param n_iter Number of iterations to obtain the distribution.
#' @param parallel If TRUE the computations will be performed in parallel.
#' @param n_cores Number of cores for the parallel computation.
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @useDynLib bamm
#' @return A data matrix of size nrow(pam) X n_iter with dispersion
#' field values.
#' @details
#'  Estimates a random distribution of the dispersion field values. To obtain
#'          random values it uses the function code{\link[bamm]{permute_pam}}
#'          at each step of the iterations. Randomization of the PAM is
#'          performed using the Babe Ruth Algorithm see Strona et al. (2014).
#'
#' @references
#' \insertRef{Soberon2015}{bamm}
#'
#' \insertRef{Strona2014}{bamm}
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' set.seed(111)
#' pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
#' dfield_rand <- bamm::null_dispersion_field_distribution(pam,n_iter=10,
#'                                                        parallel=FALSE,
#'                                                        n_cores = 2)
#' head(dfield_rand)
#' @export

null_dispersion_field_distribution <- function(pam,n_iter=10,
                                               parallel=TRUE,n_cores=2){

  if(is.data.frame(pam))
    pam <- data.matrix(pam)
  if(!is.matrix(pam))
    stop("m should be a matrix or a data.frame")
  sniter <- 1:n_iter
  if(parallel){
    plan(tweak(multisession, workers = n_cores))
  } else{
    plan(sequential)
  }

  nms <-paste0("dfrand_",sniter)
  distfield_rand <- sniter %>% furrr::future_map_dfc(function(x){
    ppam <- bamm::permute_pam(m=pam,as_sparse=TRUE)
    distfield <-bamm::pam2bioindex(pam=ppam,
                                  biodiv_index = "dispersion_field",
                                  as_sparse = FALSE)

    y <- data.frame(dfield =distfield@dispersion_field)
    return(y)
  },.progress = TRUE,.options = furrr::furrr_options(seed = NULL))
  plan(sequential)

  distfield_rand <- data.matrix(distfield_rand)
  colnames(distfield_rand) <- nms
  return(distfield_rand)
}
