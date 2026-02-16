#' csim2pam: Converts community simulation to a Presence Absence Matrix (PAM)
#' @description Converts community simulation object into a
#' Presence Absence Matrices (PAM) for a given simulation steps.
#'
#' @param community_sim An object of class \code{\link[bamm]{community_bam}}.
#' @param which_steps Steps in the simulation object to be converted into a PAM
#' @return An object of class \code{\link[bamm]{pam}}; it contains five slots.
#' 1) pams: a list of sparse matrices with Presence-Absence information (PAMs).
#' 2) which_steps: time steps corresponding to each PAM. 3) sp_names: a
#' vector of species names. 4) the grid area used in the simulation. 5) Non NA
#' cell (pixel) IDs.
#' @details For details about the object community_sim see
#' \code{\link[bamm]{community_sim}}
#' @references
#' \insertRef{SoberonOsorio}{bamm}.
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @examples
#' \donttest{
#' lagos_path <- system.file("extdata/conejos",
#'                           package = "bamm")
#' enm_path <- list.files(lagos_path,
#'                        pattern = ".tif",
#'                        full.names = TRUE)[seq(1,10)]
#' en_models <- raster::stack(enm_path)
#' ngbs_vect <- sample(1:2,replace = TRUE,
#'                     size = raster::nlayers(en_models))
#' init_coords <- read.csv(file.path(lagos_path,
#'                                   "lagos_initit.csv"))[seq(1,10),]
#' nsteps <- 10
#' sdm_comm <- bamm::community_sim(en_models = en_models,
#'                                ngbs_vect = ngbs_vect,
#'                                init_coords = init_coords,
#'                                nsteps = nsteps,
#'                                threshold = 0.1)
#'
#' pamt10 <- bamm::csim2pam(community_sim = sdm_comm ,
#'                         which_steps = 10)
#' pams <- bamm::csim2pam(community_sim = sdm_comm ,
#'                        which_steps = seq_len(10))
#' rich_pam <- bamm::pam2richness(pams,which_steps = c(1,5))
#' print(rich_pam)
#' }
#'
csim2pam <- function(community_sim, which_steps){
  if(!inherits(community_sim, "community_sim"))
    stop("Object should be of class community_sim")

  n_sps <- length(community_sim@community_sim)
  sps_names <- names(community_sim@community_sim)
  xys <- raster::xyFromCell(community_sim@community_sim[[1]]@niche_model,
                            community_sim@community_sim[[1]]@cellIDs)

  pamL <- lapply(which_steps, function(t_step){
    sim_t <- Matrix::t(community_sim@community_sim[[1]]@sdm_sim[[t_step]])
    sim_t <- Matrix::cbind2(xys,sim_t)
    for (sps in 2:n_sps) {
      sim_tm <- community_sim@community_sim[[sps]]@sdm_sim[[t_step]]
      sim_t <- Matrix::cbind2(sim_t,Matrix::t(sim_tm))
    }
    return(sim_t)
  })

  names(pamL) <- paste0("time_step_",which_steps)

  pamobj <- pam(pams =pamL,
                which_steps = which_steps,
                sp_names = sps_names,
                grid =  community_sim@community_sim[[1]]@niche_model,
                cellIDs = community_sim@community_sim[[1]]@cellIDs)


  return(pamobj)

}
