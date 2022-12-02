#' pam2richness: Converts Presence Absence Matrix (pam object) to
#' richness raster
#' @description Converts Presence Absence Matrix (pam object) to
#' richness raster
#'
#' @param pamobj An object of class pam see \code{\link[bamm]{csim2pam}}
#' @param which_steps Time steps in the pam to convert
#' @examples
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
#'                           ngbs_vect = ngbs_vect,
#'                           init_coords = init_coords,
#'                           nsteps = nsteps,
#'                           threshold = 0.1)
#'
#' pams <-bamm::csim2pam(community_sim = sdm_comm ,
#'                 which_steps = c(1:20))
#' richness_stack <- bamm::pam2richness(pams,which_steps=pams@which_steps)
#' raster::plot(richness_stack)
#' }
#' @export
#'
pam2richness <- function(pamobj,which_steps){
  if(!methods::is(pamobj,"pam")){
    stop("'pamobj' should be of class pam")
  }
  tsteps <- paste0("time_step_",which_steps)
  pams2covert <- pamobj@pams[tsteps]
  nsteps <- length(tsteps)
  pb <- utils::txtProgressBar(min = 0,
                              max = nsteps,
                              style = 3)

  richneesL <- 1:nsteps %>% purrr::map(function(x){
    grid_base <- pamobj@grid
    grid_base[pamobj@cellIDs] <- Matrix::rowSums(pams2covert[[x]])
    utils::setTxtProgressBar(pb, x)
    return(grid_base)
  })

  richness <- raster::stack(richneesL)
  names(richness) <- tsteps
  return(richness)
}
