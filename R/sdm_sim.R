#' sdm_sim: Simulate single species dispersal dynamics using the BAM framework.

#' @param set_A A setA object returned by the function
#' \code{\link[bamm]{model2sparse}}
#' @param set_M A setM object containing the adjacency matrix of the
#' study area.
#' See \code{\link[bamm]{adj_mat}}
#' @param initial_points A sparse vector returned by the function
#' \code{\link[bamm]{occs2sparse}}
#' @param nsteps Number of steps to run the simulation
#' @param stochastic_dispersal Logical. If dispersal depends on a probability of
#' visiting neighbor cells (Moore neighborhood).
#' @param disper_prop Probability of dispersal to reachable cells.
#' @param disp_prop2_suitability Logical. If probability of dispersal
#' is proportional to the suitability of reachable cells. The proportional
#' value must be declared in the parameter `disper_prop`.
#' @param progress_bar Show progress bar
#' @importFrom magrittr %>%
#' @return An object of class \code{\link[bamm]{bam}} with simulation results.
#' The simulation are stored in the sdm_sim slot (a list of sparse matrices).
#'
#' @details The model is cellular automata where the occupied area
#' of a species at time \eqn{t+1} is estimated by the multiplication of two
#' binary matrices: one matrix represents movements (M), another
#' abiotic -niche- tolerances (A)
#' (Soberon and Osorio-Olvera, 2022).
#' \deqn{\mathbf{G}_j(t+1) =\mathbf{A}_j(t)\mathbf{M}_j
#' \mathbf{G}_j(t)}
#' The equation describes a very simple process: To find the occupied patches
#' in \eqn{t+1} start with those occupied at time \eqn{t} denoted by
#' \eqn{\mathbf{G}_j(t)}, allow the individuals to disperse among
#' adjacent patches, as defined by \eqn{\mathbf{M}_j}, then remove individuals
#' from patches that are unsuitable, as defined by \eqn{\mathbf{A}_j(t)}.
#' @references
#' \insertRef{SoberonOsorio}{bamm}.
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' \donttest{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#'
#' sparse_mod <- bamm::model2sparse(model,threshold=0.05)
#' adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1)
#' occs_lep_cal <- data.frame(longitude = c(-110.08880,
#'                                          -98.89638),
#'                            latitude = c(30.43455,
#'                                         25.19919))
#'
#' occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
#'                                 occs = occs_lep_cal)
#' sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
#'                              set_M = adj_mod,
#'                              initial_points = occs_sparse,
#'                              nsteps = 10,
#'                              stochastic_dispersal = TRUE,
#'                              disp_prop2_suitability=TRUE,
#'                              disper_prop=0.5,
#'                              progress_bar=TRUE)
#'
#' sim_res <- bamm::sim2Raster(sdm_lep_cal)
#' raster::plot(sim_res)
#'
#'}
#' @export

sdm_sim <- function(set_A,set_M,initial_points,nsteps,
                    stochastic_dispersal = TRUE,
                    disp_prop2_suitability=TRUE,
                    disper_prop=0.5,progress_bar=TRUE){
  if(!inherits(set_A,"setA")){
    stop("set_A should be of class setA")
  }
  if(!inherits(set_M,"setM")){
    stop("set_M should be of class setM")
  }
  M_mat <- set_M@adj_matrix
  AMA <- set_A@sparse_model %*% M_mat # %*% set_A@sparse_model
  nsteps <- nsteps+1
  iter_vec <- 1:nsteps
  #iter_vec <- 1:50
  if(progress_bar){
    pb <- utils::txtProgressBar(min = 0,
                                max = nsteps,
                                style = 3)
  }

  g0 <- Matrix::t(initial_points)
  colnames(g0) <- set_A@cellIDs
  sdm <- list(g0)

  if(stochastic_dispersal){
    rd_adlist <- set_M@adj_list
    if(disp_prop2_suitability){
      suit_probs <- set_A@suit_values*disper_prop
      names(suit_probs) <- set_A@cellIDs
    }
    for (i in iter_vec) {
      #g0_temp <- g0
      pix_occ <- .nonzero(g0)[,2]
      nbgmat <- do.call(rbind,rd_adlist[pix_occ])
      nbgmat <- nbgmat[!duplicated(nbgmat[,2]),]
      cellIDs <- nbgmat[,2]
      n_vec <- length(cellIDs )
      if(n_vec >0){
        if(disp_prop2_suitability){
          s_probs <- suit_probs[nbgmat[,4]]
          invadible <- stats::rbinom(n = n_vec,1,
                                     prob =  s_probs )
        } else{
          invadible <- stats::rbinom(n = n_vec,1,
                                     prob = disper_prop)
        }
        M_mat[nbgmat[,3:4]] <- invadible
      }
      AMA <- set_A@sparse_model %*% M_mat #%*% set_A@sparse_model
      #M_mat <- set_M@adj_matrix
      g0 <- g0 %*% AMA
      #g0 <- g0 + g0_temp
      g0[g0>1] <- 1
      sdm[[i+1]] <- g0
      if(progress_bar){
        utils::setTxtProgressBar(pb, i)
      }

    }

  } else{
    for (i in iter_vec) {
      g0 <- g0 %*% AMA
      g0[g0>1] <- 1

      #g0 <- (g0 - min(g0))/(max(g0)-min(g0))
      sdm[[i+1]] <- g0
      if(progress_bar){
        utils::setTxtProgressBar(pb, i)
      }
    }
  }
  bam_sim <- bam(sdm_sim =sdm,
                 suit_threshold = set_A@suit_threshold,
                 niche_model=set_A@niche_model,
                 cellIDs=set_A@cellIDs,
                 sparse_model = set_A@sparse_model,
                 coordinates =set_A@coordinates,
                 adj_matrix = set_M@adj_matrix,
                 initial_points = initial_points,
                 sim_steps = nsteps-1)
  return(bam_sim)
}
#' Helper function to compute the elements in g0
#' that have no zero values.The function is taken from the
#' Ringo package
#' @param x A matrix of class "dgCMatrix"

.nonzero <- function(x){
  stopifnot(inherits(x, "dgCMatrix"))
  if (all(x@p == 0))
    return(matrix(0, nrow=0, ncol=2,
                  dimnames=list(character(0), c("row","col"))))
  res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
  colnames(res) <- c("row", "col")
  res <- res[x@x != 0, , drop = FALSE]
  return(res)
}
