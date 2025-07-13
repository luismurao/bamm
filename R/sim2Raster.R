#' sim2Raster: Convert a BAM simulation object to RasterStack
#'
#' @description Convert a BAM simulation object to RasterStack.
#' @param sdm_simul A bam object. See \code{\link[bamm]{sdm_sim}}
#' @param which_steps A numeric vector indicating the simulation steps that
#'  are going to be converted into raster layers.
#' @param progress_bar Show progress bar
#' @return A RasterStack of species' distribution at each simulation step
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @examples
#' \donttest{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)

#' sparse_mod <- bamm::model2sparse(model,threshold=0.1)
#' adj_mod <- bamm::adj_mat(sparse_mod,ngbs = 1)

#' occs_lep_cal <- data.frame(longitude = c(-115.10417,
#'                                          -104.90417),
#'                            latitude = c(29.61846,
#'                                         29.81846))

#' occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
#'                                 occs = occs_lep_cal)
#' sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
#'                             set_M = adj_mod,
#'                             initial_points = occs_sparse,
#'                             nsteps = 10)
#' sdm_lep_cal_st <- bamm::sim2Raster(sdm_simul = sdm_lep_cal,
#'                                   which_steps = seq(1,10,by=1))
#'
#' raster::plot(sdm_lep_cal_st)
#'
#' }

sim2Raster <- function(sdm_simul,which_steps = NULL,progress_bar= TRUE){

  if(!inherits(sdm_simul,"bam")){
    stop("sdm_simul should be of class bam")
  }
  if(is.null(which_steps)){
    stepsvec <- 1:sdm_simul@sim_steps
    which_steps <- stepsvec
  }
  if(max(which_steps)>sdm_simul@sim_steps)
    stop("One or more steps were not simulated by the sdm_sim function")
  if(!all(which_steps%% 1 %in% 0))
    stop("which_steps should be integers")
  else{
    stepsvec <- which_steps
  }
  grid_base <- sdm_simul@niche_model * 0
  global_vars <- c("grid_base","sdm_simul","stepsvec","progress_bar","pb")
  # Progress bar
  if (progress_bar) pb <- utils::txtProgressBar(0, max(stepsvec), style = 3)
  sim_stack <- seq_along(stepsvec) |> furrr::future_map(function(x){
    grid_base[sdm_simul@cellIDs] <- sdm_simul@sdm_sim[[stepsvec[x]]]
    if (progress_bar) utils::setTxtProgressBar(pb, x)
    return(grid_base)
  },.progress = progress_bar,
  .options = furrr::furrr_options(seed = NULL,globals = global_vars))
  sim_stack <- raster::stack(sim_stack)
  names(sim_stack) <- paste0("sim_",stepsvec)
  return(sim_stack)
}
