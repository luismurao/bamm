library(testthat)        # load testthat package
library(bamm)       # load our package
testthat::context_start_file("check-output")  # Our file is called "test-check_output.R"

# Test whether the output is a data frame
test_that("adj_mat() returns an object of class setM", {
  data("wrld_simpl", package = "maptools")
  mx <- wrld_simpl[wrld_simpl$NAME=="Mexico",]
  mx_grid <- shape2Grid(mx,0.5)
  mx_sparse <- bamm::model2sparse(mx_grid)
  adj_mx <- bamm::adj_mat(mx_sparse,ngbs=1)
  # Adjacency matrix from a niche model
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)

  sparse_mod <- bamm::model2sparse(model,threshold=0.05)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1)
  expect_s4_class(adj_mod, "setM")
})


test_that("bam_clusters() returns an object of class csd", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  model <- model > 0.7
  clusterin <- bamm::bam_clusters(model,ngbs=1,plot_model=F)
  expect_s4_class(clusterin, "csd")
})
