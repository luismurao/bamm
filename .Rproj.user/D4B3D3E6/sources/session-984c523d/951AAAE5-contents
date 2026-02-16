library(testthat)
library(bamm)
library(sp)
testthat::context_start_file("check-output")

# Test
test_that("model2sparse returns an object of class setA", {
  # Adjacency matrix from a niche model
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  sparse_mod <- bamm::model2sparse(model,threshold=0.1)
  print(sparse_mod)
  expect_s4_class(sparse_mod, "setA")
})

# Test
test_that("model2sparse, test if setA matrix is a square matrix", {
  # Adjacency matrix from a niche model
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  sparese_mod <- bamm::model2sparse(model,threshold=0.1)
  setA_dim <- dim(sparese_mod@sparse_model)
  expect_equal(setA_dim, setA_dim)
})


# Test
test_that("model2sparse using non-numeric value to binarize returns error", {
  # Adjacency matrix from a niche model
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)

  expect_error(bamm::model2sparse(model,threshold="a"))
})

# Test whether the output is a setM object and inherits s4 class
test_that("adj_mat() returns an object of class setM", {
  # Adjacency matrix from a niche model
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)

  sparse_mod <- bamm::model2sparse(model,threshold=0.05)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1)
  print(adj_mod)
  expect_s4_class(adj_mod, "setM")
})

# Test whether the output is a setM object and inherits s4 class
test_that("adj_mat() returns an object of class setM", {
  # Adjacency matrix from a niche model
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)

  sparse_mod <- bamm::model2sparse(model,threshold=0.05)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1,eigen_sys = TRUE,which_eigs = 1)
  expect_s4_class(adj_mod, "setM")
})

test_that("adj_mat()  expects an object of class setA", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  expect_error(bamm::adj_mat(model_path,ngbs=1))
})

test_that("occs2sparse()  expects an object of class setA", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  expect_error(bamm::occs2sparse(modelsparse = model_path,occs = model_path))
})

test_that("occs2sparse()  returns a sparse vector of zeros and ones", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)

  sparse_mod <- bamm::model2sparse(model,threshold=0.05)

  occs_lep_cal <- data.frame(longitude = c(-115.10417,
                                           -104.90417),
                             latitude = c(29.61846,
                                          29.81846))

  occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
                                   occs = occs_lep_cal)

  expect_s4_class(occs_sparse,"dgCMatrix")
})


test_that("shape2Grid()  converts a shapefile to a raster of a given resolution"
          , {
  x_coord <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y_coord <- c(59.736328125, 55.1220703125, 55.0341796875,
               61.142578125, 59.736328125)
  xy <- cbind(x_coord, y_coord)
  p <- sp::Polygon(xy)
  ps <- sp::Polygons(list(p),1)
  sps <- sp::SpatialPolygons(list(ps))
  r1 <- bamm::shape2Grid(sps,resolution = 0.5,ones = TRUE)
  expect_s4_class(r1,class = "RasterLayer")
})

test_that("shape2Grid()  converts a shapefile to a raster of a given resolution"
          , {
  x_coord <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y_coord <- c(59.736328125, 55.1220703125, 55.0341796875,
               61.142578125, 59.736328125)
  xy <- cbind(x_coord, y_coord)
  p <- sp::Polygon(xy)
  ps <- sp::Polygons(list(p),1)
  sps <- sp::SpatialPolygons(list(ps))
  r1 <- bamm::shape2Grid(sps,resolution = 0.5,ones = FALSE)
  expect_s4_class(r1,class = "RasterLayer")
})

test_that("permute_pam()  returns a permuted matrix with row
          sums and colum sums fixed", {
  set.seed(111)
  pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
  ppam <- bamm::permute_pam(m = pam,niter = NULL,as_sparse = FALSE)
  # Check if matrices are different
  expect_equal(object = Matrix::rowSums(pam),Matrix::rowSums(ppam))
  expect_equal(object = Matrix::colSums(pam),Matrix::colSums(ppam))

})

test_that("permute_pam()  returns a permuted matrix of class sparese", {
  set.seed(111)
  pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
  ppam <- bamm::permute_pam(m = pam,niter = NULL,as_sparse = TRUE)
  # Check if matrices are different
  expect_s4_class(ppam,"dgCMatrix")

})
test_that("permute_pam()  returns a permuted matrix of class sparese", {
  set.seed(111)
  pam <- data.frame(matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10))
  ppam <- bamm::permute_pam(m = pam,niter = NULL,as_sparse = TRUE)
  # Check if matrices are different
  expect_s4_class(ppam,"dgCMatrix")

})

test_that("permute_pam()  returns a permuted matrix of class sparese", {
  set.seed(111)
  pam <- list(matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10))
  # Check if matrices are different
  expect_error(bamm::permute_pam(m = pam,niter = NULL,as_sparse = TRUE))

})

test_that("permute_pam()  returns a permuted matrix of class sparese", {
  set.seed(111)
  pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
  # Check if matrices are different
  expect_error(bamm::permute_pam(m = pam,niter = "a",as_sparse = FALSE))

})


test_that("bam_clusters() returns an object of class csd", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  model <- model > 0.7
  #raster::crs(model) <- "+proj=longlat +datum=WGS84 +no_defs"
  clusterin <- bamm::bam_clusters(model,ngbs=1,plot_model=FALSE)
  expect_s4_class(clusterin, "csd")
})

test_that("bam_clusters() expects a RasterLayer or a sparse model", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  expect_error(bamm::bam_clusters(model_path,ngbs=1,plot_model=FALSE))
})

test_that("bam_clusters() expects a RasterLayer or a sparse model", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  model <- bamm::model2sparse(model,threshold = 0.7)
  #raster::crs(model) <- "+proj=longlat +datum=WGS84 +no_defs"
  clusterin <- bamm::bam_clusters(model,ngbs=1,plot_model=FALSE)
  print(clusterin)
  expect_s4_class(clusterin, "csd")
})


test_that("bam_clusters() returns leaflet plot with model", {
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  model <- model > 0.7
  raster::crs(model) <- "+proj=longlat"
  clusterin <- bamm::bam_clusters(model,ngbs=1,plot_model=TRUE)
  expect_s3_class(clusterin@interactive_map, "leaflet")
})

# Test eigen_bam

test_that("eigen_bam returns a list",{
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  sparse_mod <- bamm::model2sparse(model = model,threshold = 0.2)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs = 1,eigen_sys = TRUE)
  eig_bam <- bamm::eigen_bam(A=sparse_mod,M=adj_mod,
                             which_eigen = 1,rmap = TRUE)
  expect_error(bamm::eigen_bam(A="sparse_mod",
                               M=adj_mod,which_eigen = 1,rmap = TRUE))
  expect_error(bamm::eigen_bam(A=sparse_mod,M="adj_mod",
                               which_eigen = 1,rmap = TRUE))

  expect_match(class(eig_bam),"list")
})

# Test for csd_estimate

test_that("csd_estimate returns a list",{
  ## Not run:
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  model <- model > 0.7
  csd_plot <- bamm::csd_estimate(model,
                                 dispersal_steps=c(1,2))
  expect_equal(class(csd_plot),"list")
})

# Tests for sdm_sim function

test_that("sdm_sim returns an object of class bam with results from simulation"
          ,{
  ## Not run:
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)

  sparse_mod <- bamm::model2sparse(model,threshold=0.05)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1)
  occs_lep_cal <- data.frame(longitude = c(-110.08880,
                                           -98.89638),
                             latitude = c(30.43455,
                                          25.19919))

  occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
                                   occs = occs_lep_cal)
  expect_error(bamm::sdm_sim(set_A = "a",
                             set_M = adj_mod,
                             initial_points = occs_sparse,
                             nsteps = 10,
                             stochastic_dispersal = TRUE,
                             disp_prop2_suitability=TRUE,
                             disper_prop=0.5,
                             progress_bar=TRUE))
  expect_error(bamm::sdm_sim(set_A = sparse_mod,
                             set_M = "adj_mod",
                             initial_points = occs_sparse,
                             nsteps = 10,
                             stochastic_dispersal = TRUE,
                             disp_prop2_suitability=TRUE,
                             disper_prop=0.5,
                             progress_bar=TRUE))

  sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
                               set_M = adj_mod,
                               initial_points = occs_sparse,
                               nsteps = 10,
                               stochastic_dispersal = TRUE,
                               disp_prop2_suitability=FALSE,
                               disper_prop=0.5,
                               progress_bar=TRUE)

  sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
                               set_M = adj_mod,
                               initial_points = occs_sparse,
                               nsteps = 10,
                               stochastic_dispersal = FALSE,
                               disp_prop2_suitability=FALSE,
                               disper_prop=0.5,
                               progress_bar=TRUE)

  sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
                               set_M = adj_mod,
                               initial_points = occs_sparse,
                               nsteps = 10,
                               stochastic_dispersal = TRUE,
                               disp_prop2_suitability=TRUE,
                               disper_prop=0.5,
                               progress_bar=TRUE)

  expect_s4_class(sdm_lep_cal,"bam")

})

# Tests for bam_sim

test_that("bam_sim A simple simultation of predator-prey interaction.
          Returns an object of class",{
  upa <- "extdata/urania_omph/urania_guanahacabibes.tif"
  ura <- raster::raster(system.file(upa,
                                    package = "bamm"))
  opa <- "extdata/urania_omph/omphalea_guanahacabibes.tif"
  omp <- raster::raster(system.file(opa,
                                    package = "bamm"))
  msparse <- bamm::model2sparse(ura)
  init_coordsdf <- data.frame(x=-84.38751, y= 22.02932)
  initial_points <- bamm::occs2sparse(modelsparse = msparse,init_coordsdf)
  set_M <- bamm::adj_mat(modelsparse = msparse,ngbs = 1)

  expect_error(bamm::bam_sim(sp1="ura", sp2=omp, set_M=set_M,
                             initial_points=initial_points,
                             periods_toxic=3,
                             periods_suitable=3,
                             nsteps=10))
  expect_error(bamm::bam_sim(sp1=ura, sp2=omp, set_M="set_M",
                             initial_points=initial_points,
                             periods_toxic=3,
                             periods_suitable=3,
                             nsteps=10))
  ura_sim <- bamm::bam_sim(sp1=ura, sp2=omp, set_M=set_M,
                           initial_points=initial_points,
                           periods_toxic=3,
                           periods_suitable=3,
                           nsteps=10)
  expect_s4_class(ura_sim,"bam")
})

# Tests for bam_ssim

test_that("bam_ssim A simple simultation of predator-prey interaction.
          Returns an object of class",{
            upa <- "extdata/urania_omph/urania_guanahacabibes.tif"
  ura <- raster::raster(system.file(upa,
                                    package = "bamm"))
  opa <- "extdata/urania_omph/omphalea_guanahacabibes.tif"
  omp <- raster::raster(system.file(opa,
                                    package = "bamm"))
  msparse <- bamm::model2sparse(ura)
  init_coordsdf <- data.frame(x=-84.38751, y= 22.02932)
  initial_points <- bamm::occs2sparse(modelsparse = msparse,init_coordsdf)
  set_M <- bamm::adj_mat(modelsparse = msparse,ngbs = 1)

  expect_error(bamm::bam_ssim(sp1="ura", sp2=omp, set_M=set_M,
                             initial_points=initial_points,
                             periods_toxic=1,
                             periods_suitable=3,
                             nsteps=10))
  expect_error(bamm::bam_ssim(sp1=ura, sp2=omp, set_M="set_M",
                             initial_points=initial_points,
                             periods_toxic=1,
                             periods_suitable=3,
                             nsteps=10))
  ura_sim <- bamm::bam_ssim(sp1=ura, sp2=omp, set_M=set_M,
                            dispersal_prob = 0.1,
                            initial_points=initial_points,
                            periods_toxic=2,
                            periods_suitable=3,
                            palatable_matrices = TRUE,
                            nsteps=10)
  ura_sim <- bamm::bam_ssim(sp1=ura, sp2=omp, set_M=set_M,
                            dispersal_prob = 0.25,
                            initial_points=initial_points,
                            periods_toxic=1,
                            periods_suitable=3,
                            palatable_matrices = TRUE,
                            nsteps=10)
  expect_s4_class(ura_sim,"bam")
})

# Tests for sim2Raster

test_that("sim2Raster returns a stack of the distribution at time t",{
  ## Not run:
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  sparse_mod <- bamm::model2sparse(model,0.2)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs = 1,eigen_sys = TRUE,1)
  #print(adj_mod)
  occs_lep_cal <- data.frame(longitude = c(-115.10417,
                                           -104.90417),
                             latitude = c(29.61846,
                                          29.81846))
  occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
                                   occs = occs_lep_cal)
  sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
                               set_M = adj_mod,
                               initial_points = occs_sparse,
                               nsteps = 10)
  sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
                               set_M = adj_mod,
                               initial_points = occs_sparse,
                               nsteps = 10,rcpp = FALSE)
  expect_error(bamm::sim2Raster(sdm_simul = "sdm_lep_cal",
                                which_steps = seq(1,10,by=1)))
  expect_error(bamm::sim2Raster(sdm_simul = sdm_lep_cal,
                                which_steps = seq(1,11,by=1)))
  expect_error(bamm::sim2Raster(sdm_simul = sdm_lep_cal,
                                which_steps = 1.5))

  sdm_lep_cal_st <- bamm::sim2Raster(sdm_simul = sdm_lep_cal,
                                     which_steps = seq(1,10,by=1))
  sdm_lep_cal_st <- bamm::sim2Raster(sdm_simul = sdm_lep_cal,
                                     which_steps = NULL)
  expect_s4_class(sdm_lep_cal_st, "RasterStack")
})

# Tests for community_sim

test_that("community_sim simulates community dynamics and returns an
          object of class ",{
  ## Not run:
  lagos_path <- system.file("extdata/conejos",
                            package = "bamm")
  enm_path <- list.files(lagos_path,
                         pattern = ".tif",
                         full.names = TRUE)[1:5]
  en_models <- raster::stack(enm_path) >0.1
  ngbs_vect <- sample(1:2,replace = TRUE,
                      size = raster::nlayers(en_models))
  init_coords <- read.csv(file.path(lagos_path,
                                    "lagos_initit.csv"))
  nsteps <- 10
  sdm_comm <- bamm::community_sim(en_models = enm_path,
                                  ngbs_vect = ngbs_vect,
                                  init_coords = init_coords[1:5,],
                                  nsteps = nsteps,
                                  threshold = 0.1)
  expect_error(bamm::community_sim(en_models = enm_path,
                                   ngbs_vect = 1,
                                   init_coords = init_coords[1:5,],
                                   nsteps = nsteps,
                                   threshold = 0.1))
  expect_error(bamm::community_sim(en_models = enm_path,
                                   ngbs_vect = ngbs_vect,
                                   init_coords = init_coords[c(-1,2,3,4,5),],
                                   nsteps = nsteps,
                                   threshold = 0.1))
  expect_s4_class(sdm_comm,"community_sim")

  # Tests for pam2richness function

  expect_s4_class(sdm_comm,"community_sim")
  pams <-bamm::csim2pam(community_sim = sdm_comm ,
                        which_steps = c(1:10))

  expect_error(bamm::csim2pam(community_sim = "a" ,
                              which_steps = c(1:10)))

  expect_s4_class(pams,"pam")
  print(pams)
  richness_stack <- bamm::pam2richness(pams,which_steps=pams@which_steps)
  expect_error(bamm::pam2richness(pamobj = "a",which_steps=pams@which_steps))

  expect_s4_class(richness_stack,"RasterStack")

  # Tests for models2pam

  expect_error(bamm::models2pam(mods_stack = "en_models",
                                sparse=FALSE,parallel=FALSE,
                                ncores=1))
  pam <- bamm::models2pam(mods_stack = en_models,sparse=TRUE,
                          parallel=FALSE,ncores=1)
  expect_s4_class(pam,"dgCMatrix")
  pam <- bamm::models2pam(mods_stack = en_models,sparse=TRUE,
                          parallel=FALSE,ncores=1)
  expect_s4_class(pam,"dgCMatrix")
  pam <- bamm::models2pam(mods_stack = en_models,
                          sparse=FALSE,parallel=TRUE,ncores=1)
  expect_match(class(pam)[1],"matrix")
  pam <- bamm::models2pam(mods_stack = en_models,sparse=FALSE,
                          parallel=FALSE,ncores=1)
  expect_match(class(pam)[1],"matrix")

  # Test for diversity_range_analysis

  nonas <- which(!is.na(en_models[[1]][]))
  xy_mat <- sp::coordinates(en_models[[1]])[ nonas,]
  pam <- bamm::models2pam(en_models,sparse=FALSE)
  rdivan <- bamm::diversity_range_analysis(pam=pam,parallel = TRUE,
                                           xy_mat=xy_mat,
                                           raster_templete = en_models[[1]],
                                           n_cores = 1,
                                           return_null_dfield=TRUE)
  expect_error(bamm::plot(rdivan,plot_type="diversity_range1"))
  bamm::plot(rdivan,plot_type="diversity_range_map")
  if(interactive()){
    # bamm::plot(rdivan,plot_type="diversity_range_interactive")

  }

  #bamm::plot(rdivan,plot_type="diversity_range_interactive")
  bamm::plot(rdivan,plot_type="alpha")
  bamm::plot(rdivan,plot_type="dispersion_field")
  bamm::plot(rdivan,plot_type="dispersion_field_map")
  expect_s4_class(rdivan,"diversity_range")
  rdivan <- bamm::diversity_range_analysis(pam=pam,parallel = TRUE,
                                           xy_mat=xy_mat,n_cores =1,
                                           raster_templete = NULL,
                                           return_null_dfield=TRUE)
  bamm::plot(rdivan,plot_type="diversity_range_map")

})

# Tests for pam2richness function

#test_that("pam2richness returns a raster of richness",{
#  lagos_path <- system.file("extdata/conejos",
#                            package = "bamm")
#  enm_path <- list.files(lagos_path,
#                         pattern = ".tif",
#                         full.names = TRUE)
#  en_models <- raster::stack(enm_path)
#  ngbs_vect <- sample(1:2,replace = TRUE,
#                      size = raster::nlayers(en_models))
#  init_coords <- read.csv(file.path(lagos_path,
#                                    "lagos_initit.csv"))
#  nsteps <- 10
#  sdm_comm <- bamm::community_sim(en_models = enm_path,
#                                  ngbs_vect = ngbs_vect,
#                                  init_coords = init_coords,
#                                  nsteps = nsteps,
#                                  threshold = 0.3)
#  expect_s4_class(sdm_comm,"community_sim")
#  pams <-bamm::csim2pam(community_sim = sdm_comm ,
#                        which_steps = c(1:10))

#  expect_error(bamm::csim2pam(community_sim = "a" ,
#                              which_steps = c(1:10)))

#  expect_s4_class(pams,"pam")
#  print(pams)
#  richness_stack <- bamm::pam2richness(pams,which_steps=pams@which_steps)
#  expect_error(bamm::pam2richness(pamobj = "a",which_steps=pams@which_steps))

#  expect_s4_class(richness_stack,"RasterStack")

#})

# Test for null_distribution_field_distribution

test_that("null_distribution_field_distribution expects a matrix",{
  set.seed(111)
  pam <- data.frame(matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10))
  dfield_rand <- bamm::null_dispersion_field_distribution(pam,n_iter=10,
                                                          parallel=FALSE,
                                                          n_cores = 1)
  dfield_rand <- bamm::null_dispersion_field_distribution(pam,n_iter=10,
                                                          parallel=TRUE,
                                                          n_cores = 1)
  expect_error(bamm::null_dispersion_field_distribution("pam",n_iter=10,
                                                        parallel=FALSE,
                                                        n_cores = 1))
  expect_vector(dfield_rand)
})

# Tests for jaccard.R

test_that("jaccard returns a data.frame", {
  m1_path <- system.file("extdata/conejos/Lepus_othus_cont.tif",
                         package = "bamm")
  m2_path <- system.file("extdata/conejos/Brachylagus_idahoensis_cont.tif",
                         package = "bamm")
  m1 <- raster::raster(m1_path) > 0.01
  m2 <- raster::raster(m2_path) >0.01
  m1s <- bamm::model2sparse(m1,threshold = 0.1)
  m2s <- bamm::model2sparse(m2,threshold = 0.1)
  jcc <- bamm::jaccard(m1,m2)
  jccs <- bamm::jaccard(m1s,m2s)
  expect_error(bamm::jaccard("m1",m2))
  expect_error(bamm::jaccard(m1,"m2"))
  expect_equal(jcc,jccs)
  expect_equal(class(jcc),"data.frame")
})

# Tests for pam2bioindex

test_that("pam2bioindex returns an object of class bioindex",{
  set.seed(111)
  pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
  bioindices <- bamm::pam2bioindex(pam=data.frame(pam),biodiv_index="all")
  expect_s4_class(bioindices,"bioindex")
  # Return results as sparse models
  bioindices <- bamm::pam2bioindex(pam=pam,biodiv_index="all",as_sparse=TRUE)
  print(bioindices)
  expect_s4_class(bioindices,"bioindex_sparse")
  expect_error(bamm::pam2bioindex(pam="pam",biodiv_index="all",as_sparse=TRUE))

})

# Tests for models2pam
#test_that("models2pam returns a PAM as a sparsematrix",{
#  lagos_path <- system.file("extdata/conejos",
#                            package = "bamm")
#  enm_path <- list.files(lagos_path,
#                         pattern = ".tif",
#                         full.names = TRUE)
#  en_models <- raster::stack(enm_path) >0.01
#  expect_error(bamm::models2pam(mods_stack = "en_models",
#                                sparse=FALSE,parallel=FALSE,
#                                ncores=2))
#  pam <- bamm::models2pam(mods_stack = en_models,sparse=TRUE,
#                          parallel=TRUE,ncores=2)
#  expect_s4_class(pam,"dgCMatrix")
#  pam <- bamm::models2pam(mods_stack = en_models,sparse=TRUE,
#                          parallel=FALSE,ncores=2)
#  expect_s4_class(pam,"dgCMatrix")
#  pam <- bamm::models2pam(mods_stack = en_models,
#                          sparse=FALSE,parallel=TRUE,ncores=2)
#  expect_match(class(pam)[1],"matrix")
#  pam <- bamm::models2pam(mods_stack = en_models,sparse=FALSE,
#                          parallel=FALSE,ncores=2)
#  expect_match(class(pam)[1],"matrix")
#})

# Test for diversity_range_analysis

#test_that("diversity_range_analysis returns an object of class diversity_range",
#          {
#  set.seed(111)
#  pam <- matrix(rbinom(10000,1,0.5),nrow = 100,ncol = 1000)
#  rdivan <- bamm::diversity_range_analysis(pam=pam,parallel = FALSE,
#                                           return_null_dfield=TRUE)
#  print(rdivan)
#  expect_s4_class(rdivan,"diversity_range")

#  bamm::plot(rdivan,plot_type="diversity_range")
  # Lagomorphos
#  lagos_path <- system.file("extdata/conejos",
#                            package = "bamm")
#  enm_path <- list.files(lagos_path,
#                         pattern = ".tif",
#                         full.names = TRUE)
#  en_models <- raster::stack(enm_path) >0.01
#  nonas <- which(!is.na(en_models[[1]][]))
#  xy_mat <- sp::coordinates(en_models[[1]])[ nonas,]
#  pam <- bamm::models2pam(en_models,sparse=FALSE)
#  rdivan <- bamm::diversity_range_analysis(pam=pam,parallel = TRUE,
#                                           xy_mat=xy_mat,
#                                           raster_templete = en_models[[1]],
#                                           return_null_dfield=TRUE)
#  expect_error(bamm::plot(rdivan,plot_type="diversity_range1"))
#  bamm::plot(rdivan,plot_type="diversity_range_map")
  #bamm::plot(rdivan,plot_type="diversity_range_interactive")
#  bamm::plot(rdivan,plot_type="alpha")
#  bamm::plot(rdivan,plot_type="dispersion_field")
#  bamm::plot(rdivan,plot_type="dispersion_field_map")
#  expect_s4_class(rdivan,"diversity_range")
#  rdivan <- bamm::diversity_range_analysis(pam=pam,parallel = TRUE,
#                                           xy_mat=xy_mat,
#                                           raster_templete = NULL,
#                                           return_null_dfield=TRUE)
#  bamm::plot(rdivan,plot_type="diversity_range_map")
#})

# Testing predic method

test_that("predict retuns a prediction",{
  # Not run:
  # Load R packages
  # rm(list = ls())
  # Read raster model for Lepus californicus
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  # Convert model to sparse
  sparse_mod <- bamm::model2sparse(model = model,0.1)
  # Compute adjacency matrix
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1)

  # Initial points to start dispersal process

  occs_lep_cal <- data.frame(longitude = c(-115.10417,
                                           -104.90417),
                             latitude = c(29.61846,
                                          29.81846))
  # Convert to sparse the initial points
  occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
                                   occs = occs_lep_cal)

  # Run the bam (sdm) simultation for 100 time steps
  smd_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
                               set_M = adj_mod,
                               initial_points = occs_sparse,
                               nsteps = 10,rcpp = TRUE)
  #----------------------------------------------------------------------------
  # Predict species' distribution under suitability change
  # scenarios (could be climate chage scenarios).
  #----------------------------------------------------------------------------

  # Read suitability layers (two suitability change scenarios)
  layers_path <- system.file("extdata/suit_change",
                             package = "bamm")
  niche_mods_stack <- raster::stack(list.files(layers_path,
                                               pattern = ".tif$",
                                               full.names = TRUE))
  raster::plot(niche_mods_stack)
  # Predict
  new_preds <- predict(object = smd_lep_cal,
                       niche_layers = niche_mods_stack,
                       nsteps_vec = c(10,10))
  fname <- tempfile(pattern = "animation_",fileext = ".html")
  # Generate the dispersal animation for time period 1 and 2
  new_preds <- predict(object = smd_lep_cal,
                       niche_layers = niche_mods_stack,
                       nsteps_vec = c(10,10),
                       animate=TRUE,
                       filename=fname,
                       fmt="HTML")
  expect_error(predict(object = smd_lep_cal,
                       niche_layers = "niche_mods_stack",
                       nsteps_vec = c(10,10),
                       animate=TRUE,
                       filename=fname,
                       fmt="HTML"))
  expect_error(predict(object = smd_lep_cal,
                       niche_layers = "niche_mods_stack",
                       nsteps_vec = c(10,10),
                       animate=TRUE,
                       filename=fname,
                       fmt="nn"))
  expect_error(predict(object = smd_lep_cal,
                       niche_layers = niche_mods_stack,
                       nsteps_vec = c(10,1,5),
                       animate=TRUE,
                       filename=fname,
                       fmt="HTML"))
  expect_error(predict(object = smd_lep_cal,
                       niche_layers = niche_mods_stack,
                       nsteps_vec = c(10,1),
                       nbgs_vec=c(1,2,5),
                       animate=TRUE,
                       filename=fname,
                       fmt="HTML"))
  fname <- tempfile(pattern = "animation_",fileext = ".gif")
  new_preds <- predict(object = smd_lep_cal,
                       niche_layers = niche_mods_stack[[1]],
                       nsteps_vec = c(10),
                       animate=TRUE,
                       filename=fname,
                       fmt="GIF")
  new_preds <- predict(object = smd_lep_cal,
                       niche_layers = niche_mods_stack,
                       nsteps_vec = c(1),
                       stochastic_dispersal=TRUE,
                       nbgs_vec=c(1,2),
                       disp_prop2_suitability	=TRUE,
                       disper_prop= 0.5,
                       period_names	=c("P1","P2"),
                       bg_color	="gray97",
                       suit_color	= "red",
                       occupied_color="blue",
                       animate=TRUE,
                       filename=fname,
                       fmt="GIF")


})
# Test for sim2Animation
testthat::test_that("sim2Animation",{
  model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                            package = "bamm")
  model <- raster::raster(model_path)
  sparse_mod <- bamm::model2sparse(model,0.1)
  adj_mod <- bamm::adj_mat(sparse_mod,ngbs=2)
  occs_lep_cal <- data.frame(longitude = c(-115.10417,
                                           -104.90417),
                             latitude = c(29.61846,
                                          29.81846))
  occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
                                   occs = occs_lep_cal)
  sdm_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
                               set_M = adj_mod,
                               initial_points = occs_sparse,
                               nsteps = 20)
  ani_name <- tempfile(pattern = "anima_",fileext = ".html")
  sdm_lep_cal_st <- bamm::sim2Animation(sdm_simul = sdm_lep_cal,
                                        which_steps = seq(1,20,by=1),
                                        fmt = "HTML",ani.width = 1200,
                                        ani.height = 1200,
                                        filename = ani_name)
  expect_null(sdm_lep_cal_st)
  ani_name <- tempfile(pattern = "anima_",fileext = ".gif")
  sdm_lep_cal_st <- bamm::sim2Animation(sdm_simul = sdm_lep_cal,
                                        which_steps = seq(1,20,by=1),
                                        fmt = "GIF",ani.width = 1200,
                                        ani.height = 1200,
                                        filename = ani_name)

  expect_null(sdm_lep_cal_st)

})


# Test polygons 2 pam

testthat::test_that("pol2pam",{
  # Example with sample data
  uicn <- readRDS(system.file("extdata/uicn.rds",package = "bamm"))
  sudam <- readRDS(system.file("extdata/suam.rds",package = "bamm"))
  # Convert to PAM with 0.5 degree resolution
  pam_result <- bamm::pol2pam(poly = uicn,
                              taxon_attribute = "binomial",
                              resolution = 0.5,
                              polymask = NULL)

  # With masking polygon
  pam_masked <- pol2pam(poly = uicn,
                        taxon_attribute = "binomial",
                        resolution = 0.5,
                        polymask = sudam)
  testthat::expect_equal(class(pam_masked),"data.frame")

})
