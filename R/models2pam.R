#' models2pam: Converts binary rasters to a PAM
#' @description Function to convert binary raster models to a
#' Presence Absences Matrix.
#' @param mods_stack A raster stack containing binary models of each
#' species in the community.
#' @param return_coords Logical. If TRUE the pam will be returned with
#' coordinates in the first two columns.
#' @param sparse Logical. If TRUE the PAM will be returned as a sparse matrix.
#' @param parallel Logical. If TRUE computations will be done in parallel
#' @param ncores Integer. Number of cores to run the parallel process.
#' @return A presence-absence matrix (PAM).
#' @details For more information about PAM see Soberon and Cavner (2015).
#' @references
#' \insertRef{Soberon2015}{bamm}.
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @examples
#' \donttest{
#' lagos_path <- system.file("extdata/conejos",
#'                           package = "bamm")
#' enm_path <- list.files(lagos_path,
#'                        pattern = ".tif",
#'                        full.names = TRUE)[1:10]
#' en_models <- raster::stack(enm_path) >0.01
#' pam <- bamm::models2pam(en_models,
#'                         return_coords=TRUE,
#'                         sparse=FALSE,
#'                         parallel=FALSE,ncores=1)
#' head(pam)
#' }
models2pam <- function(mods_stack,return_coords=FALSE,sparse=TRUE,
                       parallel=FALSE,ncores=2){
  cmod <-class(mods_stack)
  if(!cmod %in% c("RasterStack","RasterBrick")){
    stop("'mods_stack' should be of class 'RasterStack' or 'RasterBrick'")
  }
  else{
    mods_stack <- raster::stack(mods_stack)
    nsps <- raster::nlayers(mods_stack)
    rvals <- mods_stack[[1]][]
    cellIDs <- which(!is.na(rvals))
    if(parallel){
      oplan <- plan(multisession,workers=ncores)
      on.exit(plan(oplan), add = TRUE)
    }

    if(!sparse){
      pam0 <- seq_len(nsps) |> furrr::future_map_dfc(function(x){
        m1 <- mods_stack[[x]]
        m2 <- m1[]*1
        m3 <- m2[cellIDs]
        df1 <- data.frame(m3)
        df1 <- stats::setNames(df1,names(m1))
        return(df1)
      },.progress = TRUE,
      .options = furrr::furrr_options(seed = NULL,
                                      globals = c("mods_stack","cellIDs"),
                                      packages = "raster"))
      pam0 <- data.matrix(pam0)
      if(return_coords){
        xys <- sp::coordinates(mods_stack[[1]])[cellIDs,]
        pam0 <- data.frame(xys,pam0)
      }
      return(pam0)
    } else{
      pamL <- seq_len(nsps) |> furrr::future_map(function(x){
        m1 <- mods_stack[[x]]
        m2 <- m1[]*1.0
        m3 <- m2[cellIDs]
        unosIDs <- which(m3==1)
        msparse0 <- Matrix::sparseVector(x = rep(100,length(unosIDs)),
                                         i = unosIDs,
                                         length=length(cellIDs))

        return(msparse0)
      },.progress = TRUE,
      .options = furrr::furrr_options(seed = NULL,
                                      globals = c("mods_stack","cellIDs"),
                                      packages = "raster"))

      pamL <- lapply(pamL, as, "sparseMatrix")
      pam0 <- pamL[[1]]
      for(i in 2:length(pamL)){
        pam0 <- Matrix::cbind2(pam0,pamL[[i]])
      }

      #pam0 <- do.call(cbind2, pamL)
      colnames(pam0) <- names(mods_stack)
      if(return_coords){
        xys <- sp::coordinates(mods_stack[[1]])[cellIDs,]
        xys <- Matrix::Matrix(xys,sparse = TRUE)
        pam0 <- Matrix::cbind2(xys,pam0)
      }
      return(pam0)
    }
  }
}

