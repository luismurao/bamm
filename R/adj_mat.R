#' adj_mat: Function to compute the adjacency matrix of an area.
#'
#' @description Creates an adjacency matrix of an area of interest.
#' This could be the accessible area (M) of a species or any geographic
#' region of interest.
#' @param modelsparse A setA object returned by the function
#' \code{\link[bamm]{model2sparse}}.
#' @param ngbs Numeric. Number of neighbors (see details).
#' @param eigen_sys Logical. If TRUE the eigen analyses of the
#' adjacency matrix will be returned.
#' @param which_eigs Numeric. Which eigen value and eigen vector will
#' be returned.
#' @return Returns an object of class \code{\link[bamm]{setM}} with 7 slots.
#' The first contains the adjacency matrix. A n x n sparse matrix (n=number of
#' non-NA cells of the niche model) where connected cells are represented by 1.
#' The second slot has the adjacency list. It is a list of matrices with four
#' columns (FromRasCell -from cell ID of the raster-, -to cell ID of the
#' raster-, -from non-NA cell-, -to non-NA cell-). Other slots contain
#' information about initial coordinates where dispersal occurs
#' (initial_points), number of cells used to define the neighborhood (ngbs),
#' non-NA coordinates (coordinates), and a matrix of eigen vectors (eigen_vec).
#' @export
#' @details The model is a raster object of the area where the dispersal
#' process will occur.
#'
#' The function creates an adjacency matrix where cells are considered connected based on:
#' \itemize{
#'   \item The specified neighborhood size (`ngbs`)
#'   \item The spatial resolution of the input raster
#'   \item Only non-NA cells in the original model
#' }
#'
#' When `eigen_sys = TRUE`, the function performs spectral decomposition using \code{\link[RSpectra]{eigs}},
#' which is particularly useful for:
#' \itemize{
#'   \item Analyzing network connectivity properties
#'   \item Identifying clusters or communities in the landscape
#'   \item Modeling dispersal processes using spectral graph theory
#' }
#'
#' The number of neighbors depends on the dispersal abilities of the species
#' and the spatial resolution of the niche model; for example, a species's
#' with big dispersal abilities will move throughout more than 1 km^2 per day,
#' so the idea is to give an approximate number of moving neighbors (pixels)
#' per unit of time.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[bamm]{model2sparse}} for creating the input `setA` object
#'   \item \code{\link[bamm]{setM}} for the output class structure
#'   \item \code{\link[RSpectra]{eigs}} for the eigen decomposition implementation
#' }
#' For more information about see adjacency matrices in the context of
#' the theory of area of distribution (Soberon and Osorio-Olvera, 2022).
#'
#' @references
#' \insertRef{SoberonOsorio}{bamm}.
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' x_coord <- c(-106.5699, -111.3737,-113.9332,
#'              -110.8913, -106.4262, -106.5699)
#' y_coord <- c(16.62661, 17.72373, 19.87618,
#'              22.50763, 21.37728, 16.62661)
#' xy <- cbind(x_coord, y_coord)
#' p <- sp::Polygon(xy)
#' ps <- sp::Polygons(list(p),1)
#' sps <- sp::SpatialPolygons(list(ps))
#' mx_grid <- bamm::shape2Grid(sps,resolution = 0.25,ones = TRUE)
#' mx_sparse <- bamm::model2sparse(model=mx_grid, threshold = 0.1)
#' adj_mx <- bamm::adj_mat(modelsparse=mx_sparse,
#'                         ngbs=1,eigen_sys=TRUE,which_eigs=1)
#' print(adj_mx)
#' mx_grid_eigen <- mx_grid
#' mx_grid_eigen[mx_sparse@cellIDs] <- adj_mx@eigen_vec
#' raster::plot(mx_grid_eigen)
#'


adj_mat <- function(modelsparse,ngbs=1,eigen_sys=FALSE,which_eigs=1){


  if(!inherits(modelsparse,"setA")){
    stop("modelsparse should be of class setA")
  }

  nbase <- 2*ngbs+1
  ngMat <- base::matrix(rep(1,nbase*nbase),
                        ncol =nbase,byrow = TRUE)
  ngMat[ngbs+1,ngbs+1] <- 0

  no_na <- modelsparse@cellIDs

  r_ad <- raster::adjacent(x = modelsparse@niche_model,
                           cells = no_na,directions = ngMat,
                           target = no_na)

  m_ad1 <- Matrix::sparseMatrix( i=match(r_ad[,1],no_na),
                                 j=match(r_ad[,2],no_na),
                                 x=1.0 )
  id_nona <- seq_along(no_na)
  newff <- as.factor(c(r_ad[,1],r_ad[,2]))
  newff2 <- newff
  connected_cells <- unique(c(r_ad[,1],r_ad[,2]))
  connected_ids <- id_nona[which(no_na %in% connected_cells)]
  levels(newff2) <- connected_ids
  newnu <- as.numeric(as.character(newff2))
  idc <- seq_len(nrow(r_ad))
  from <- as.numeric(as.character(newff[idc]))
  to <- as.numeric(as.character(newff[-idc]))
  from_nu <- newnu[idc]
  to_nu <- newnu[-idc]
  big_vec <- c(from,to, from_nu,to_nu)
  r_ad_b <- matrix(big_vec,ncol = 4,byrow = FALSE)
  colnames(r_ad_b) <- c("FromRasCell","ToRasCell",
                        "FromNonNaCell","ToNonNaCell")
  rd_adlist <- split.data.frame(r_ad_b, r_ad_b[,3])

  g_set0 <- setM(adj_matrix = m_ad1,
                 adj_list = rd_adlist,
                 ngbs =ngbs,
                 coordinates = modelsparse@coordinates)

  if(eigen_sys){
    eigSys <- RSpectra::eigs(A = m_ad1,k=which_eigs)
    g_set0 <- setM(adj_matrix = m_ad1,
                   coordinates = modelsparse@coordinates,
                   adj_list = rd_adlist,
                   ngbs =ngbs,
                   eigen_vec = eigSys$vectors,
                   eigen_val = eigSys$values)
  }



  return(g_set0)
}

