#' csim2pam: Converts community simulation to a Presence Absence Matrix (PAM)
#' @description Converts community simulation object into a
#' Presence Absence Matrices (PAM) for a given simulation steps.
#'
#' @param community_sim An object of class community_bam.
#' @param which_steps Steps in the simulation object to be converted into a PAM
#' @export
#' @examples
#'
#' \dontrun{
#' lagos_path <- system.file("extdata/conejos",
#'                           package = "bamm")
#' enm_path <- list.files(lagos_path,
#'                        pattern = ".tif",
#'                     full.names = TRUE)
#' en_models <- raster::stack(enm_path)
#' ngbs_vect <- sample(1:2,replace = T,
#'                     size = raster::nlayers(en_models))
#' init_coords <- read.csv(file.path(lagos_path,
#'                                   "lagos_initit.csv"))
#' nsteps <- 30
#' sdm_comm <- bamm::community_sim(en_models = enm_path,
#'                                ngbs_vect = ngbs_vect,
#'                                init_coords = init_coords,
#'                                nsteps = nsteps,
#'                                threshold = 0.1)
#'
#' pamt10 <- bamm::csim2pam(community_sim = sdm_comm ,
#'                        which_steps = 10)
#' pams <- bamm::csim2pam(community_sim = sdm_comm ,
#'                       which_steps = c(1:10))
#' }
#'
csim2pam <- function(community_sim, which_steps){
  if(!inherits(community_sim, "community_sim"))
    stop("Object should be of class community_sim")

  n_sps <- length(community_sim@community_sim)
  sps_names <- names(community_sim@community_sim)
  pamL <- lapply(which_steps, function(t_step){
    sim_t <- Matrix::t(community_sim@community_sim[[1]]@sdm_sim[[t_step]])
    for (sps in 2:n_sps) {
      sim_tm <- community_sim@community_sim[[sps]]@sdm_sim[[t_step]]
      sim_t <- Matrix::cbind2(sim_t,Matrix::t(sim_tm))
      #dimnames(sim_t) <- list(sps_names)
    }
    #rows=paste0("s",1:nrow(sim_t))
    #columns <- paste0(sps_names,1:ncol(sim_t))
    #dimnames(sim_t ) = list(rows,columns)
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
