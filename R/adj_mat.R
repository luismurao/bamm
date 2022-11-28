#' adj_mat: Function to compute the adjacency matrix of an area.
#'
#' @description Creates an adjacency matrix of an area of interest. This could be the accessible area (M) of a species or any geographic region of interest.
#' @param modelsparse A setA object returned by the function \code{\link[bamm]{model2sparse}}.
#' @param ngbs Numeric. Number of neighbors (see details).
#' @param eigen_sys Logical. If TRUE the eigen analsys of the adjancency matrix will be returned.
#' @param which_eigs Numeric. Which eigen value and eigen vector will be returned.
#' @return Returns an adjacency matrix of class sparseMatrix of n+2 x n columns (n number of the non-NA cells of grid_base) with the coordinates of the non-NA cells of grid_base.
#' @export
#' @details The grid_base raster object is the area where the dispersal process will occur.
#' The number of neighbors depends on the dispersal abilities of the species and the spatial resolution
#' of the grid_base; for example, a species's with big dispersal abilities will move throughout more
#' than 1 km^2 per day, so the idea is to give an approximate number of moving neighbors (pixels) per
#' unit of time.
#'
#' @examples
#' \dontrun{
#' data("wrld_simpl", package = "maptools")
#' mx <- wrld_simpl[wrld_simpl$NAME=="Mexico",]
#' mx_grid <- shape2Grid(mx,0.5)
#' mx_sparse <- bamm::model2sparse(mx_grid)
#' adj_mx <- bamm::adj_mat(mx_sparse,ngbs=1)
#' # Adjacency matrix from a niche model
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#'
#' sparse_mod <- bamm::model2sparse(model,threshold=0.05)
#' adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1)
#' }


adj_mat <- function(modelsparse,ngbs=1,eigen_sys=FALSE,which_eigs=1){


  if(!inherits(modelsparse,"setA")){
    stop("modelsparse should be of class setA")
  }

  nbase <- 2*ngbs+1
  ngMat <- base::matrix(rep(1,nbase*nbase),
                        ncol =nbase,byrow = T )
  ngMat[ngbs+1,ngbs+1] <- 0

  no_na <- modelsparse@cellIDs

  r_ad <- raster::adjacent(x = modelsparse@niche_model,
                           cells = no_na,directions = ngMat,
                           target = no_na)

  m_ad1 <- Matrix::sparseMatrix( i=match(r_ad[,1],no_na),
                                 j=match(r_ad[,2],no_na),
                                 x=1.0 )
  id_nona <- 1:length(no_na)
  newff <- as.factor(c(r_ad[,1],r_ad[,2]))
  newff2 <- newff
  connected_cells <- as.numeric(as.character(levels(newff)))
  connected_ids <- id_nona[which(no_na %in% connected_cells)]
  levels(newff2) <- connected_ids
  newnu <- as.numeric(as.character(newff2))
  idc <- 1:nrow(r_ad)
  from <- as.numeric(as.character(newff[idc]))
  to <- as.numeric(as.character(newff[-idc]))
  from_nu <- newnu[idc]
  to_nu <- newnu[-idc]
  big_vec <- c(from,to, from_nu,to_nu)
  r_ad_b <- matrix(big_vec,ncol = 4,byrow = F)
  rd_adlist <- split.data.frame(r_ad_b, r_ad_b[,3])

  g_set0 <- setM(adj_matrix = m_ad1,
                 adj_list = rd_adlist,
                 ngbs =ngbs,
                 coordinates = modelsparse@coordinates)

  if(eigen_sys){
    eigSys <- RSpectra::eigs(A = m_ad1,k=which_eigs)
    g_set0 <- setM(adj_matrix = m_ad1,
                   coordinates = modelsparse@coordinates,
                   eigen_vec = eigSys$vectors,
                   eigen_val = eigSys$values)
  }



  return(g_set0)
}



# Modifed: A slightly modification of the function
# .adjacentUD (in raster package) from Robert Hijmans
# Licence GPL v3



.adjacentBAM <- function(x, cells, ngb) {

  directions <- sum(ngb==1, na.rm=TRUE)

  rs <- raster::res(x)
  target <- cells
  rn <- raster::raster(ngb)
  center <- which(raster::values(rn)==0)
  rc <- raster::rowFromCell(rn, center)
  cc <- raster::colFromCell(rn, center)
  xngb <- yngb <- ngb
  xngb[] <- rep(1:ncol(ngb), each=nrow(ngb)) - cc
  yngb[] <- rep(nrow(ngb):1, ncol(ngb)) - (nrow(ngb)-rc+1)
  ngb[ngb != 1] <- NA
  xngb <- stats::na.omit(as.vector( xngb * rs[1] * ngb))
  yngb <- stats::na.omit(as.vector( yngb * rs[2] * ngb))

  xy <- raster::xyFromCell(x, cells)
  X <- unlist(lapply(as.vector(xy[,1,drop=FALSE]), function(z) z + xngb ))
  Y <- unlist(lapply(as.vector(xy[,2,drop=FALSE]), function(z) z + yngb ))

  d <- data.frame(X,Y)


  cell <- rep(cells, each=directions)

  d <- stats::na.omit(cbind(cell, raster::cellFromXY(x, d)))

  attr(d, 'na.action') <- NULL
  colnames(d) <- c('fromCell', 'toCell')
  d <- d[d[,2] %in% target, ]

  d <- d[order(d[,1], d[,2]),]

  return(d)
}


