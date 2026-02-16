#' sdm_sim: Simulate single species dispersal dynamics using the BAM framework.

#' @param set_A A setA object returned by the function
#' \code{\link[bamm]{model2sparse}}
#' @param set_M A setM object containing the adjacency matrix of the
#' study area.
#' See \code{\link[bamm]{adj_mat}}
#' @param sp_interaction_model A RasterLayer representing a suitability model
#' of a species positive interaction. See details
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
#' @param rcpp Logical. Use native C++ code to run the model.
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
#'
#' The stochastic model uses a suitability values to model dispersal probabilities.
#' These suitability values can be either obtained from the continues model stored in
#' the set_A object or from the sp_interaction_model (RasterLayer).
#' If the parameter rcpp is set to TRUE the model will be run using native C++
#' code through Rcpp and RcppArmadillo packages.
#'
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
#'                              nsteps = 100,
#'                              stochastic_dispersal = TRUE,
#'                              disp_prop2_suitability=FALSE,
#'                              disper_prop=0.5,
#'                              progress_bar=TRUE)
#'
#' sim_res <- bamm::sim2Raster(sdm_lep_cal)
#' raster::plot(sim_res)
#'
#'}
#' @export

sdm_sim <- function(set_A, set_M,sp_interaction_model=NULL ,initial_points, nsteps,
                    stochastic_dispersal = TRUE,
                    disp_prop2_suitability = TRUE,
                    disper_prop = 0.5, progress_bar = TRUE,rcpp=TRUE) {

  # Validate inputs
  if (!inherits(set_A, "setA")) stop("set_A should be of class setA")
  if (!inherits(set_M, "setM")) stop("set_M should be of class setM")
  suit_probs <- set_A@suit_values
  nsteps <- nsteps + 1

  if(!is.null(sp_interaction_model)){
    if(inherits(sp_interaction_model,"RasterLayer")){
      test_mod <- raster::stack(set_A@niche_model,sp_interaction_model)
      suit_probs <- sp_interaction_model[set_A@cellIDs]
      rm(test_mod)
    } else{
      stop("sp_interaction_model should be of class RasterLayer")
    }
  }

  if(rcpp){
    # Extract components in R
    A_sparse <- as(set_A@sparse_model, "dgCMatrix")
    M_adj <- as(set_M@adj_matrix, "dgCMatrix")
    adj_list <- set_M@adj_list

    # Call C++ function with extracted components
    sdm <- sdm_sim_rcpp(
      A = A_sparse,
      M_orig = M_adj,
      g0_input = initial_points,
      suit_values = suit_probs,
      adj_list = adj_list,
      nsteps = nsteps,
      stochastic_dispersal = stochastic_dispersal,
      disp_prop2_suitability = disp_prop2_suitability,
      disper_prop = disper_prop,
      progress_bar = progress_bar
    )

  } else{
    # BAM objects
    M_orig <- set_M@adj_matrix
    AMA <- set_A@sparse_model %*% M_orig %*% set_A@sparse_model
    g0 <- initial_points
    rownames(g0) <- set_A@cellIDs
    sdm <- list(g0)

    # Progress bar
    if (progress_bar) pb <- utils::txtProgressBar(0, nsteps, style = 3)

    if (stochastic_dispersal) {
      rd_adlist <- set_M@adj_list

      # Precompute suitability probabilities if needed
      if (disp_prop2_suitability) {
        suit_probs <-suit_probs * disper_prop
        names(suit_probs) <- set_A@cellIDs
      }

      for (i in seq_len(nsteps)) {
        # Reset to original matrix (critical fix!)
        M_mat <- M_orig

        # Get occupied cells
        occ_indices <- .nonzero(g0)[, 2]

        if (length(occ_indices) > 0) {
          # Efficiently process neighbors
          all_edges <- lapply(occ_indices, function(idx) {
            if (is.null(rd_adlist[[idx]])) return(NULL)
            rd_adlist[[idx]]
          })

          if (!is.null(all_edges)) {
            nbgmat <- do.call(rbind, all_edges)
            n_vec <- nrow(nbgmat)

            # Calculate dispersal probabilities
            if (disp_prop2_suitability) {
              s_probs <- suit_probs[nbgmat[, 4]]
              invadible <- stats::rbinom(n_vec, 1, s_probs)
            } else {
              invadible <- stats::rbinom(n_vec, 1, disper_prop)
            }

            # Update adjacency matrix
            M_mat[nbgmat[, 3:4]] <- invadible
            M_mat[nbgmat[, 4:3]] <- invadible
            #M_mat <- Matrix::forceSymmetric(M_mat)  # Maintain symmetry
          }
        }

        # Propagate dispersal
        AMA <- set_A@sparse_model %*% M_mat  %*% set_A@sparse_model
        g0 <-  AMA %*% g0
        g0 <-  g0 + sdm[[i]]
        g0@x[g0@x > 1] <- 1  # Efficient thresholding

        sdm[[i + 1]] <- g0
        if (progress_bar) utils::setTxtProgressBar(pb, i)
      }

    } else {
      # Deterministic dispersal
      for (i in seq_len(nsteps)) {
        g0 <-  AMA %*% g0
        g0@x[g0@x > 1] <- 1  # Efficient thresholding
        sdm[[i + 1]] <- g0
        if (progress_bar) utils::setTxtProgressBar(pb, i)
      }
    }

    if (progress_bar) close(pb)
  }



  # Return bam object
  bam(
    sdm_sim = sdm,
    suit_threshold = set_A@suit_threshold,
    niche_model = set_A@niche_model,
    cellIDs = set_A@cellIDs,
    sparse_model = set_A@sparse_model,
    coordinates = set_A@coordinates,
    adj_matrix = set_M@adj_matrix,
    adj_list = set_M@adj_list,
    initial_points = initial_points,
    sim_steps = nsteps - 1
  )
}
#' Helper function to compute the elements in g0
#' that have no zero values.The function is taken from the
#' Ringo package
#' @param x A matrix of class "dgCMatrix"
#' @noRd

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
