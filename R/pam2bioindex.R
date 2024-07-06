#' pam2bioindex: PAM to biodiversity index
#' @description pam2bioindex estimates various biodiversity indices for a
#' certain PAM.
#' @param pam A Presence-Absence-Matrix of matrix class or sparse matrix.
#' @param biodiv_index Possible values are alpha, omega,
#'  wbeta (Whittaker’s multiplicative beta index),
#'  laBeta (Lande’s additive beta index)
#'  dispersion_field, all.
#' @param as_sparse Return indices as sparse objects
#' @importFrom Rdpack reprompt
#' @return An object of class \code{\link[bamm]{bioindex}} with three slots
#' each represents a matrix of diversity indices: alpha, omega, and
#' dispersion field, richness_field.
#' @details The biodiversity indices can be found in Soberón and Cavner (2015).
#' @references
#' \insertRef{Soberon2015}{bamm}
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @examples
#' set.seed(111)
#' pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
#' bioindices <- bamm::pam2bioindex(pam=pam,biodiv_index="all")
#' # Return results as sparse models
#' bioindices <- bamm::pam2bioindex(pam=pam,biodiv_index="all",as_sparse=TRUE)
#' bioindices@alpha
#' bioindices@omega
#' bioindices@dispersion_field
#' @useDynLib bamm

pam2bioindex <- function(pam,biodiv_index="dispersion_field",as_sparse=FALSE) {

  if(is.data.frame(pam)){
    pam <- as.matrix(pam)
    pam <- Matrix::Matrix(pam, sparse = TRUE)
  }
  if(is.matrix(pam)){
    pam <- Matrix::Matrix(pam, sparse = TRUE)
  }
  if(!methods::is(pam,"dgCMatrix")){
    stop('pam should be of class "matrix" or "dgCMatrix"')
  }
  # N the number of sites
  N <- nrow(pam)
  # S the number of species
  S <- ncol(pam)

  if(as_sparse)
    bioindices <- methods::new(Class = "bioindex_sparse")
  if(!as_sparse)
    bioindices <- methods::new(Class="bioindex")

  if(any(c("alpha","richness_field","all") %in% biodiv_index)){
    ones_sps <- rep(1,S)
    alpha <- pam %*% ones_sps
    bioindices@alpha <- if(as_sparse) alpha else as(alpha,"matrix")
  }
  if(any(c("omega","dispersion_field",
           "wbeta","labeta","lebeta","all") %in% biodiv_index)){
    ones_sites <- rep(1,N)
    pamT <- Matrix::t(pam)
    om   <- pamT  %*% ones_sites
    trace_om <- sum(om)
    bioindices@omega <- if(as_sparse) om else as(om,"matrix")
  }
  if(any(c("wbeta","all") %in% biodiv_index)){
    wbeta <- S*N/trace_om
    bioindices@wBeta <- wbeta

  }
  if(any(c("labeta","all") %in% biodiv_index)){
    labeta <- S*(1-(trace_om/(S*N)))
    bioindices@laBeta <- labeta
  }

  if(any(c("dispersion_field","lebeta","all") %in% biodiv_index)){
    om_st <- om/N
    fist <- pam%*%om_st
    bioindices@dispersion_field <- if(as_sparse) fist else as(fist,"matrix")
  }
  if(any(c("richness_field", "all") %in% biodiv_index)){
    ones_sps <- rep(1,S)
    pamT <- Matrix::t(pam)
    richf <- pamT %*% alpha
    bioindices@richness_field <- if(as_sparse) richf else as(richf,"matrix")
  }
  if(any(c("lebeta","all") %in% biodiv_index)){
    ones_sites <- rep(1,N)
    #fist <- pam%*%om_st
    #lebeta <- trace_om - Matrix::t(fist) %*%  ones_sites
    lebeta <- ( S/wbeta)* Matrix::t(fist) %*%  ones_sites
    bioindices@leBeta <- as.numeric(lebeta)
  }
  if(any(c("nestedness", "all") %in% biodiv_index)){
    ones_sites <- rep(1,N)
    fist_ <- pam%*%om
    nss <- (Matrix::t(fist_) %*% ones_sites) - (N*S/wbeta)
    bioindices@nestedness <- as.numeric(nss)
  }

  return(bioindices)
}
