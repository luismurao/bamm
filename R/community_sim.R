#' community_bam: Community \code{bamm}
#'
#' @description Estimate community dynamics using the \code{bamm} framework
#' @param en_models A stack or directory with the ecological niche models for
#' each species in the community.
#' @param ngbs_vect A vector
#' @param init_coords A data.frame with 3 columns: sp_name, x and
#' y; x is the longitude and y is the latitude of initial dispersal points
#' @param nsteps Number of iteration steps for the simulation.
#' @param parallel TRUE if the simulation should be run in parallel.
#' Default FALSE
#' @param ... Parameters of the functions of the \code{bamm} class
#' @export
#' @examples
#' \dontrun{
#' lagos_path <- system.file("extdata/conejos",
#'                           package = "bamm")
#' enm_path <- list.files(lagos_path,
#'                        pattern = ".tif",
#'                        full.names = TRUE)

#' en_models <- raster::stack(enm_path)

#' ngbs_vect <- sample(1:2,replace = TRUE,
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
#'}

community_sim <- function(en_models,ngbs_vect,init_coords,nsteps,
                          parallel, ...){

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
  sdmsimL <- 1:n_models %>% purrr::map(function(x){

    #sparse_mod <- bamm::model2sparse(model = models[[x]],
    #                                ...)
    sparse_mod <- bamm::model2sparse(model = models[[x]],
                                    threshold = 0.1)
    # Compute adjacency matrix
    adj_mod <- bamm::adj_mat(sparse_mod,ngbs=ngbs_vect[[x]])

    # Initial points to start dispersal process
    init_coords <- init_coordsL[[x]][,2:3]
    occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
                                    occs = init_coords)

    sdmsimul <- bamm::sdm_sim(set_A = sparse_mod,
                             set_M = adj_mod,
                             initial_points = occs_sparse,
                             nsteps = nsteps)
    utils::setTxtProgressBar(pb, x)

    return(sdmsimul)

  })
  names(sdmsimL) <- names(init_coordsL)

  res1 <- community_bam(community_sim = sdmsimL)

  return(res1)
}
