#' null_dispersion_field_distribution: Null distribution of the
#' dispersion field
#' @description null_dispersion_field_distribution estimates a
#' random distribution of the dispersion field values.
#' @param pam A Presence-Absence-Matrix of matrix class or sparse matrix.
#' @param n_iter Number of iterations to obtain the distribution.
#' @param parallel If TRUE the computations will be performed in parallel.
#' @param n_cores Number of cores for the parallel computation.
#' @param randal Randomization algorithm applied to the PAM.
#' Possible choices "curveball", "fastball", and "indep_swap".
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @useDynLib bamm
#' @return A data matrix of size nrow(pam) X n_iter with dispersion
#' field values.
#' @details
#'  Estimates a random distribution of the dispersion field values. To obtain
#'          random values it uses the function \code{\link[bamm]{permute_pam}}
#'          at each step of the iterations. Randomization of the PAM can be
#'          performed using the "fastball" (Godard and Neal, 2022) and the
#'          "curveball" (Strona et al., 2014), and  and the independent
#'          swap (Kembel et al. 2010) algorithms.
#'          The implementation of the "fastball" in C++ is provided
#'          in \url{https://github.com/zpneal/fastball/blob/main/fastball.cpp}
#'
#' @references
#' \insertRef{Soberon2015}{bamm}
#'
#' \insertRef{Strona2014}{bamm}
#'
#' \insertRef{Gordard2022}{bamm}
#'
#' \insertRef{Kembel2010}{bamm}
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' set.seed(111)
#' pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
#' dfield_rand <- bamm::null_dispersion_field_distribution(pam,n_iter=10,
#'                                                        parallel=FALSE,
#'                                                        randal="indep_swap",
#'                                                        n_cores = 2)
#' head(dfield_rand)
#' @export

null_dispersion_field_distribution <- function(pam,n_iter=10,randal="fastball",
                                               parallel=TRUE,n_cores=2){
  ral <- match.arg(arg = randal,
                   choices = c("indep_swap","curveball","fastball"))

  if(is.data.frame(pam))
    pam <- data.matrix(pam)
  if(!methods::is(pam,"matrix") & !is.numeric(pam[1,1])){
    stop("pam object should be a binary matrix")
  }
  sniter <- 1:n_iter
  if(parallel){
    oplan <- plan(tweak(multisession, workers = n_cores))
  } else{
    oplan <- plan(sequential)
  }
  on.exit(plan(oplan), add = TRUE)
  nms <-paste0("dfrand_",sniter)
  distfield_rand <- sniter %>% furrr::future_map_dfc(function(x){
    ppam <- bamm::permute_pam(m=pam,as_sparse=TRUE,randal = randal)
    distfield <-bamm::pam2bioindex(pam=ppam,
                                   biodiv_index = "dispersion_field",
                                   as_sparse = FALSE)

    y <- data.frame(dfield =distfield@dispersion_field)
    return(y)
  },.progress = TRUE,.options = furrr::furrr_options(seed = NULL))


  distfield_rand <- data.matrix(distfield_rand)
  colnames(distfield_rand) <- nms
  if(randal == "fastball"){
    message("This function uses the code from:\n",
            "https://github.com/zpneal/fastball/blob/main/fastball.cpp\n",
            "Please when using the function for publications cite:\n
             Godard K, Neal ZP (2022). 'fastball: a fast algorithm to randomly
             sample bipartite graphs with fixed degree sequences'.
             Journal of Complex Networks, 10(6), cnac049.")
  }
  return(distfield_rand)
}
