#' community_bam: Community \code{bamm}
#'
#' @description Estimate community dynamics using the \code{bamm} framework
#' @param en_models A stack or directory with the ecological niche models for
#' each species in the community.
#' @param ngbs_vect A vector containing the number of neighbors for each
#' adjacency matrix of each species in the community
#' see \code{\link[bamm]{adj_mat}}.
#' @param init_coords A data.frame with 3 columns: sp_name, x and
#' y; x is the longitude and y is the latitude of initial dispersal points
#' @param nsteps Number of iteration steps for the simulation.
#' @param threshold_vec A vector of threshold values used to bnarize niche
#' models.
#' @param stochastic_dispersal Logical. If dispersal depends on a probability of
#' visiting neighbor cells (Moore neighborhood).
#' @param disper_prop Probability of dispersal to reachable cells.
#' @param disp_prop2_suitability Logical. If probability of dispersal
#' is proportional to the suitability of reachable cells. The proportional
#' value must be declared in the parameter `disper_prop`.
#' @return An object of class community_sim. The object contains simulation
#' results for each species in the community.
#' @details Each element in community_sim is an object of class. For more
#' details about the simulation see \code{\link[bamm]{sdm_sim}}.
#' \code{\link[bamm]{bam}}.
#' @references
#' \insertRef{SoberonOsorio}{bamm}.
#' @author Luis Osorio-Olvera & Jorge Soberon
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
#' nsteps <- 12

#' sdm_comm <- bamm::community_sim(en_models = en_models,
#'                                 ngbs_vect = ngbs_vect,
#'                                 init_coords = init_coords,
#'                                 nsteps = nsteps)
#'
#' com_pam <- bamm::csim2pam(sdm_comm,which_steps = seq(1,nsteps))
#' rich_pam <- pam2richness(com_pam,which_steps = c(1,5,10))
#' raster::plot(rich_pam)
#'}

community_sim <- function(en_models,
                          ngbs_vect,
                          init_coords,
                          nsteps,
                          threshold_vec = NULL,
                          stochastic_dispersal = FALSE,
                          disp_prop2_suitability=TRUE,
                          disper_prop=0.5){

  models <- raster::stack(en_models)
  n_models <- raster::nlayers(models)
  n_ngbs_vect <- length(ngbs_vect)
  if(n_models != n_ngbs_vect){
    stop(paste0("Each model should have a dispersal hypothesis,",
                "in other words ngbs_vect should have the same",
                "number of elements as the number of niche models"))
  }
  init_coordsL <- split(init_coords,init_coords[[1]])

  if(length(init_coordsL) != n_ngbs_vect){
    stop(paste0("Each model should have initial points,",
                "in other words init_coords should have the same",
                "number of species as the number of niche models"))
  }

  pb <- utils::txtProgressBar(min = 0,
                              max = n_models,
                              style = 3)
  if(is.null(threshold_vec) || length(threshold_vec) != n_models){
    threshold_vec <- rep(0.1,n_models)
  }
  sdmsimL <- seq_len(n_models) %>% furrr::future_map(function(x){

    sparse_mod <- bamm::model2sparse(model = models[[x]],
                                     threshold = threshold_vec[x])
    #sparse_mod <- bamm::model2sparse(model = models[[x]],
    #                                threshold = 0.1)
    # Compute adjacency matrix
    adj_mod <- bamm::adj_mat(sparse_mod,ngbs=ngbs_vect[[x]])

    # Initial points to start dispersal process
    init_coords <- init_coordsL[[x]][,2:3]
    occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
                                     occs = init_coords)

    sdmsimul <- bamm::sdm_sim(set_A = sparse_mod,
                              set_M = adj_mod,
                              initial_points = occs_sparse,
                              nsteps = nsteps,
                              stochastic_dispersal = stochastic_dispersal,
                              disp_prop2_suitability=disp_prop2_suitability,
                              disper_prop=disper_prop)
    utils::setTxtProgressBar(pb, x)

    return(sdmsimul)

  },.progress = TRUE,.options = furrr::furrr_options(seed = NULL))
  names(sdmsimL) <- names(init_coordsL)

  res1 <- community_bam(community_sim = sdmsimL)

  return(res1)
}

