## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bamm)

## -----------------------------------------------------------------------------
model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                          package = "bamm")
model <- raster::raster(model_path) 
# binary model
model_bin <- model > 0.7
sparse_mod <- bamm::model2sparse(model = model_bin)


## -----------------------------------------------------------------------------
sparse_mod

## -----------------------------------------------------------------------------
 # Adjacency matrix from a niche model


adj_mod <- adj_mat(sparse_mod,ngbs=1,eigen_sys = T)
adj_mod

## ----fig.cap="Fig. 1 Connectivity map of `adj_mod` object. Here green color means more connected."----
model_eig <- model
model_eig[sparse_mod@cellIDs] <- abs(adj_mod@eigen_vec)
raster::plot(model_eig)

## -----------------------------------------------------------------------------
clustersin <- bamm::bam_clusters(model=sparse_mod,
                                ngbs=1,plot_model=FALSE)

## -----------------------------------------------------------------------------
clustersin

## ----fig.cap="Figure 2. An interative map showing the geographic clusters for a species that can travel two steps per unit time"----
clustersin@interactive_map

## ----fig.cap="Figure 3. An raster map showing the geographic clusters for a species that can travel two steps per unit time"----
raster::plot(clustersin@raster_map)

## ----fig.width=4,fig.height=4,fig.cap="Figure 4. Connectivity Sutitability Dispersal plot (CSD plot).\nThe mean number of connected cells (MNCC) is showed in the legend."----
csd_plot <- bamm::csd_estimate(sparse_mod,
                              dispersal_steps=c(2,4,8))
plot(csd_plot$csd$`dispersal_step:_4`@raster_map,
     main="CSD clusters")


