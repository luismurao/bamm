#' pam2richness: Converts Presence Absence Matrix (pam object) to
#' richness raster
#' @description Converts Presence Absence Matrix (pam object) to
#' richness raster
#'
#' @param pamobj An object of class pam see \code{\link[bamm]{csim2pam}}
#' @param which_steps Time steps in the pam to convert
#' @return A RasterStack richness for each simulation step
#'
#' @author Luis Osorio-Olvera & Jorge Soberón.
#' @examples
#' \donttest{
#' lagos_path <- system.file("extdata/conejos",
#'                           package = "bamm")
#' enm_path <- list.files(lagos_path,
#'                        pattern = ".tif",
#'                     full.names = TRUE)[seq(1,10)]
#' en_models <- raster::stack(enm_path)
#' ngbs_vect <- sample(2,replace = TRUE,
#'                     size = raster::nlayers(en_models))
#' init_coords <- read.csv(file.path(lagos_path,
#'                                   "lagos_initit.csv"))[seq(1,10),]
#' nsteps <- 10
#' sdm_comm <- bamm::community_sim(en_models = en_models,
#'                                 ngbs_vect = ngbs_vect,
#'                                 init_coords = init_coords,
#'                                 nsteps = nsteps,
#'                                 threshold = 0.1)
#'
#' pams <-bamm::csim2pam(community_sim = sdm_comm ,
#'                       which_steps = seq_len(nsteps))
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
