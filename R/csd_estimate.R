#' csd_estimate: Estimate the connectivity suitability and dispersal plot
#' @description csd_plot gives an estimate of the number of geographic clusters
#' given a set of dispersal hypothesis and a suitability raster
#'
#' @param model A raster model or a setA object representing the
#' suitability model
#' @param dispersal_steps A numeric vector with elements representing
#' the dispersal hypothesis to test.
#' @return A list of length three. The first element contains the Connectivity-
#' Suitability-Diagram information estimated for each element in the vector
#' of dispersal_steps. The second is tbl_df object with a summary of the number
#' of cluster of each dispersal step and the mean number of connected clusters.
#' The last element is base plot showing the information cointained in
#' the tbl_df object.
#' @details For more information about the Connectivity-Suitability-Diagram
#' see \code{\link[bamm]{bam_clusters}}
#' @references
#' \insertRef{SoberonOsorio}{bamm}.
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @examples
#' \donttest{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#' model <- model > 0.7
#' csd_plot <- bamm::csd_estimate(model,
#'                          dispersal_steps=c(2,4,8))
#' csd_plot$plot
#' }
#' @export

csd_estimate <- function(model,dispersal_steps=c(2,4,8,16,32,64)){

  if(methods::is(model,"RasterLayer")){
    model <- bamm::model2sparse(model)
  }
  dispersal_steps <- sort(dispersal_steps)
  ds <- seq_along(dispersal_steps)
  pb <- utils::txtProgressBar(min = 0, max = max(ds), style = 3)
  #testworks <- TRUE
  csd <- ds %>% purrr::map(function(x){
    #if(testworks){
    bclust <- try({
      r <- bamm::bam_clusters(model,ngbs = dispersal_steps[x])
       },silent = TRUE)

      #if(!methods::is(bclust,"csd")){
      #  warning("There is not enough memory to calculate",
      #          "the adjacency matrix for dispersal step =",
      #          dispersal_steps[x],"\n returning ",
      #          paste(dispersal_steps[1:(x-1)],collapse = " "))
      #  testworks <- FALSE
      #}
      #}
      utils::setTxtProgressBar(pb, x)
      return(bclust)
  })



  ndisp <- length(csd)
  d_all <- 1:ndisp %>% purrr::map_df(function(x){
    if(methods::is(csd[[x]],"csd")){
      d1 <- data.frame(csd[[x]]@connections,
                       d=dispersal_steps[x])
      return(d1)
    }
  })

  #d_clust <- d_all %>%
  #  dplyr::group_by(d) %>%
  #  dplyr::summarise(n_clusters=max(clusterID))
  d <- NULL; clusterID <- NULL; nclust <-NULL
  d_clust <- d_all  %>% dplyr::group_by(d,clusterID) %>%
    dplyr::summarise(nclust=dplyr::n()) %>% dplyr::group_by(d) %>%
    dplyr::summarise(Clusters=max(clusterID),
                     mean_area=mean(nclust))

  s1 <- as.factor(d_clust$mean_area)
  nlevs <- length(unique(d_clust$mean_area))
  levels(s1) <- seq(1,2,length.out = nlevs)
  #data.frame(f1=s1,
  #           f2=d_clust$mean_area)
  sizes <- as.numeric(as.character(s1))

  graphics::plot(d_clust$d,d_clust$Clusters,frame=TRUE,
                 type="l",pch=19,lwd=2,cex=2,
                 xlab="Neighbors",ylab="Clusters",
                 cex.lab=1.25, cex.axis=1.25)
  graphics::points(d_clust$d,d_clust$Clusters,pch=19,
                   cex=sizes*1.5,col="brown")


  #sizes <- ( d_clust$mean_area - min(d_clust$mean_area))/
  # (max(d_clust$mean_area)-min(d_clust$mean_area))
  #sizes[1] <- sizes[2]/2
  #sizes <- sizes*5


  graphics::legend("topright",
                   legend = paste(round(d_clust$mean_area)),
                   col = "brown",
                   pch = 19,
                   bty = "n",
                   pt.cex = sizes,
                   cex = 1,
                   title="MNCC",
                   text.col = "black",
                   horiz = FALSE,
                   inset = c(0.1, 0.1))

  p <- grDevices::recordPlot()
  #graphics::plot.new() ## clean up device
  #
  names(csd) <- paste0("dispersal_step:_",dispersal_steps)
  csd <- list(csd=csd,plot_data=d_clust,plot=p)


  return(csd)

}





