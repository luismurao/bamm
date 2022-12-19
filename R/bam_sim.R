#' bam_sim: Simulate dispersal dynamics using the set B of the BAM framework.

#' @param sp1 Niche model of the focal species (the one that disperses).
#' @param sp2 Niche model of the species with whom sp1 interacts
#' (currently no dispersal dynamics for this species).
#' @param set_M A setM object containing the adjacency matrix for sp1.
#' See \code{\link[bamm]{adj_mat}}
#' @param periods_toxic  Time periods that sps2 takes to develop defense
#' mechanisms (i.e. toxic).
#' @param periods_suitable This is the time that sp2 takes to become non-toxic
#' @param initial_points A sparse vector returned by the function
#' \code{\link[bamm]{occs2sparse}}
#' @param nsteps Number of steps to run the simulation
#' @param progress_bar Show progress bar
#' @importFrom magrittr %>%
#' @return An object of class bam. The object contains 12 slots of information
#' (see details) from which simulation results are stored in sdm_sim object,
#' a list of sparse matrices with results of each simulation step.
#' @details The returned object inherits from \code{\link[bamm]{setA}},
#' \code{\link[bamm]{setM}} classes. Details about the dynamic model
#' can be found in Soberon and Osorio-Olvera (2022). The model is cellular
#' automata where the occupied area of a species at time \eqn{t+1} is
#' estimated by the multiplication of three
#' binary matrices: one matrix represents movements (M), another
#' abiotic -niche- tolerances (A), and a third, biotic interactions (B)
#' (Soberon and Osorio-Olvera, 2022).
#' \deqn{\mathbf{G}_j(t+1) =\mathbf{B}_j(t)\mathbf{A}_j(t)\mathbf{M}_j
#' \mathbf{G}_j(t)}
#' @references
#' \insertRef{SoberonOsorio}{bamm}.
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' # Compute dispersal dynamics of Urania boisduvalii as a function of
#' # palatable Omphalea
#' urap <- system.file("extdata/urania_omph/urania_guanahacabibes.tif",
#'                                   package = "bamm")
#' ura <- raster::raster(urap)
#' ompp <- system.file("extdata/urania_omph/omphalea_guanahacabibes.tif",
#'                                   package = "bamm")
#' omp <- raster::raster(ompp)
#' msparse <- bamm::model2sparse(ura)
#' init_coordsdf <- data.frame(x=-84.38751, y= 22.02932)
#' initial_points <- bamm::occs2sparse(modelsparse = msparse,init_coordsdf)
#' set_M <- bamm::adj_mat(modelsparse = msparse,ngbs = 1)
#' ura_sim <- bamm::bam_sim(sp1=ura, sp2=omp, set_M=set_M,
#'                          initial_points=initial_points,
#'                          periods_toxic=5,
#'                          periods_suitable=1,
#'                          nsteps=40)
#' ura_omp <- bamm::sim2Raster(ura_sim)
#' raster::plot(ura_omp[[c(1,5,10,15,20,30,35,40)]])
#' \donttest{
#' if(requireNamespace("animation")){
#' # Animation example
#' anp <-tempfile(pattern = "simulation_results_",fileext = ".gif")
#' new_sim <- bamm::sim2Animation(sdm_simul = ura_sim,
#'                               which_steps = seq_len(ura_sim@sim_steps),
#'                               fmt = "GIF",
#'                               filename = anp)
#'}
#'}
#' @export
#'
#'

bam_sim <- function(sp1,sp2,set_M,initial_points,
                     periods_toxic,periods_suitable,
                     nsteps,progress_bar=TRUE){
  if(!(methods::is(sp1,"RasterLayer") && methods::is(sp2,"RasterLayer"))){
    stop("sp1 and sp2 should be of raster class")
  }
  if(!inherits(set_M,"setM")){
    stop("set_M should be of class setM")
  }
  sp1_sp2 <- sp1*sp2
  A <- bamm::model2sparse(sp1_sp2)
  bin_model <- A@sparse_model
  Matrix::diag(set_M@adj_matrix) <- 1
  AMA <- A@sparse_model %*% set_M@adj_matrix %*%  A@sparse_model
  AMA2 <- AMA
  g0 <- initial_points
  nonzerosL <- list()
  time_mat <- matrix(numeric(length(A@cellIDs)),ncol  = 1)
  time_counter_off <- Matrix::Matrix(data = time_mat,sparse=TRUE)
  time_counter_on <- Matrix::Matrix(data = time_mat,sparse=TRUE)

  sdm <- list(g0)
  if(progress_bar){
    pb <- utils::txtProgressBar(min = 0,
                                max = nsteps,
                                style = 3)
  }

  for(i in 1:nsteps){
    pix_occ <- .nonzero(g0)[,1]
    g0 <- AMA%*%g0
    g0[g0>1] <- 1
    time_counter_off[pix_occ, ] <- time_counter_off[pix_occ, ] + 1
    to_off_vec <- pix_occ[which( time_counter_off[pix_occ, ] >= periods_toxic)]
    time_onID <- .nonzero(time_counter_on)[,1]

    if(length(to_off_vec)>0L){
      g0[to_off_vec] <- 0
      Matrix::diag(bin_model)[to_off_vec] <- 0
      AMA <- bin_model %*% set_M@adj_matrix
      time_counter_off[to_off_vec, ] <- 0
      time_onID <- c(time_onID,to_off_vec)
    }
    time_counter_on[time_onID, ] <- time_counter_on[time_onID, ]  + 1
    to_on_vec <- time_onID[which(time_counter_on[time_onID,]>=periods_suitable)]
    if(length(to_on_vec)>0L){
      Matrix::diag(bin_model)[to_on_vec] <- 1
      AMA <- bin_model %*% set_M@adj_matrix
      time_counter_on[to_on_vec,] <- 0
    }
    sdm[[i+1]] <- g0
    if(progress_bar){
      utils::setTxtProgressBar(pb, i)
    }

  }

  bamsim <- bam(sdm_sim =sdm,
                niche_model=sp1_sp2,
                cellIDs=A@cellIDs,
                sparse_model = A@sparse_model,
                coordinates =A@coordinates,
                adj_matrix = set_M@adj_matrix,
                initial_points = initial_points,
                sim_steps = nsteps)
  return(bamsim)

}


