#' bam_clusters: Function to estimate the connectivity of suitable areas
#'
#' @description Function to estimate the connectivity of suitable areas given
#' an adjacency matrix.
#' @param model A niche model in raster format
#' @param ngbs Numeric. Number of neighbors (see details).
#' @param plot_model Logical. Indicates whether to plot the model in the
#' cluster map.
#' @details
#' The grid_base raster object is the area where the dispersal process will
#' occur.
#' The number of neighbors depends on the dispersal abilities of the species
#' and the spatial resolution of the grid_base;
#' for example, a species's with big dispersal abilities will move throughout
#' more than 1 km^2 per day, so the idea is to give an approximate number of
#' moving neighbors (pixels) per unit of time.
#' @importFrom methods as
#' @importFrom raster as.factor
#' @return A list with a data.frame of the coordinates of each cluster and a
#' leaflet map.
#'
#' @examples

#'
#' \dontrun{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#' model <- model > 0.7
#' clusterin <- bamm::bam_clusters(model,ngbs=1,plot_model=TRUE)
#' clusterin@interactive_map
#' }
#' @export
#'
#'
#'
#'
bam_clusters <- function(model,ngbs=1,plot_model=FALSE){

  if(methods::is(model,"RasterLayer")){
    msparse <- bamm::model2sparse(model)
  }
  else if(methods::is(model,"setA")){
    msparse <- model
  }
  else{
    stop("model should be a RasterLayer or setA class object")
  }
  #if(class(model) != "setA")
  #  stop("model should be of raster class or setA class")


  adj_msp <- bamm::adj_mat(msparse,ngbs)
  ama <-   msparse@sparse_model %*% adj_msp@adj_matrix %*%  msparse@sparse_model

  ids_p <- which(Matrix::rowSums(ama)>0)
  ama <- ama[ids_p ,ids_p]
  colnames(ama) <- paste0(ids_p)
  rownames(ama) <-  paste0(ids_p)

  ama <- as(ama, "TsparseMatrix")
  rgs <- ama@i + 1L
  cls <- ama@j + 1L
  vals <- ama@x
  my_df <- data.frame(
    rows = factor(ama@Dimnames[[1]][rgs],
                  levels=unique(ama@Dimnames[[1]])),
    cols = factor(ama@Dimnames[[2]][cls],
                  levels=unique(ama@Dimnames[[2]])),
    value = vals)
  net <- igraph::graph.data.frame(my_df, directed = FALSE)
  cl <- igraph::clusters(net)

  to_find <- seq_along(cl$csize)[cl$csize > 1]

  clusterDF <- to_find %>%
    purrr::map_df(function(x){
      cellID <- as.numeric(igraph::V(net)$name[cl$membership %in% x])
      dcl <- data.frame(cluster=x,
                        cluster_size = cl$csize[x],
                        cellID= cellID)
      return(dcl)
    })

  paleta <- grDevices::rainbow(length(to_find)*500)



  # Mapa de los clusters


  raster_clust <- msparse@niche_model

  df_clust1 <- data.frame(msparse@coordinates[clusterDF$cellID,],
                          clusterID = clusterDF$cluster,
                          cluster_size=clusterDF$cluster_size)

  raster_clust[msparse@cellIDs[clusterDF$cellID]] <- clusterDF$cluster

  cluster_map <- paste("cluster: ",as.character(clusterDF$cluster),
                       "<br>","c_size: ",clusterDF$cluster_size,
                       #"<br>","p_ID: ",clusterDF$cellID,
                       "<br>","Longitude: ",round(df_clust1$x,4),
                       "<br>","Latitude: ", round(df_clust1$y,4))
  ids_cols <- floor(seq(1,length(paleta),
                  length.out = max(clusterDF$cluster)))

  cluster_color <- paleta[ids_cols]
  cluster_color <- sample(cluster_color)

  cols <- cluster_color[clusterDF$cluster]

  m <- leaflet::leaflet(df_clust1) %>% leaflet::addTiles()

  if(plot_model){
    #nw <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    #m2 <- raster::projectRaster(model, crs=nw)
    #m2 <- leaflet::projectRasterForLeaflet(model,method = "bilinear")

    mod <- round(msparse@niche_model)
    mod <- raster::as.factor(mod)
    #raster::crs(mod)
    m <- m %>%
      leaflet::addRasterImage( mod,
                               colors = c("gray100","blue"),
                               opacity = 0.5)
  }

  m <- m %>% leaflet::addCircleMarkers(lng = ~x,
                                       lat = ~y,
                                       popup = cluster_map,radius = 0.1,
                                       color = cols,opacity = 1)

  csd_res <- csd(connections = df_clust1,
                 interactive_map = m,
                 raster_map = raster_clust)
  #attributes(csd_res)$interactive_map <- m

  return(csd_res)

}
