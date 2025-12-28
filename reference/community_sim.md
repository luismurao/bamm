# community_bam: Community `bamm`

Estimate community dynamics using the `bamm` framework

## Usage

``` r
community_sim(
  en_models,
  ngbs_vect,
  init_coords,
  nsteps,
  threshold_vec = NULL,
  stochastic_dispersal = FALSE,
  disp_prop2_suitability = TRUE,
  disper_prop = 0.5,
  rcpp = TRUE
)
```

## Arguments

- en_models:

  A stack or directory with the ecological niche models for each species
  in the community.

- ngbs_vect:

  A vector containing the number of neighbors for each adjacency matrix
  of each species in the community see [`adj_mat`](adj_mat.md).

- init_coords:

  A data.frame with 3 columns: sp_name, x and y; x is the longitude and
  y is the latitude of initial dispersal points

- nsteps:

  Number of iteration steps for the simulation.

- threshold_vec:

  A vector of threshold values used to bnarize niche models.

- stochastic_dispersal:

  Logical. If dispersal depends on a probability of visiting neighbor
  cells (Moore neighborhood).

- disp_prop2_suitability:

  Logical. If probability of dispersal is proportional to the
  suitability of reachable cells. The proportional

- disper_prop:

  Probability of dispersal to reachable cells.

- rcpp:

  Logical. Use native C++ code to run the model. value must be declared
  in the parameter \`disper_prop\`.

## Value

An object of class community_sim. The object contains simulation results
for each species in the community.

## Details

Each element in community_sim is an object of class. For more details
about the simulation see [`sdm_sim`](sdm_sim.md). [`bam`](bam-class.md).

## References

Soberón J, Osorio-Olvera L (2023). “A dynamic theory of the area of
distribution.” *Journal of Biogeography6*, **50**, 1037-1048.
[doi:10.1111/jbi.14587](https://doi.org/10.1111/jbi.14587) ,
<https://onlinelibrary.wiley.com/doi/abs/10.1111/jbi.14587>. .

## Author

Luis Osorio-Olvera & Jorge Soberon

## Examples

``` r
# \donttest{
lagos_path <- system.file("extdata/conejos",
                          package = "bamm")
enm_path <- list.files(lagos_path,
                       pattern = ".tif",
                       full.names = TRUE)[seq(1,10)]
en_models <- raster::stack(enm_path)
ngbs_vect <- sample(1:2,replace = TRUE,
                    size = raster::nlayers(en_models))
init_coords <- read.csv(file.path(lagos_path,
                                  "lagos_initit.csv"))[seq(1,10),]
nsteps <- 12
sdm_comm <- bamm::community_sim(en_models = en_models,
                                ngbs_vect = ngbs_vect,
                                init_coords = init_coords,
                                nsteps = nsteps)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%

com_pam <- bamm::csim2pam(sdm_comm,which_steps = seq(1,nsteps))
#>   |                                                                              |                                                                      |   0%  |                                                                              |======                                                                |   8%  |                                                                              |============                                                          |  17%  |                                                                              |==================                                                    |  25%  |                                                                              |=======================                                               |  33%  |                                                                              |=============================                                         |  42%  |                                                                              |===================================                                   |  50%  |                                                                              |=========================================                             |  58%  |                                                                              |===============================================                       |  67%  |                                                                              |====================================================                  |  75%  |                                                                              |==========================================================            |  83%  |                                                                              |================================================================      |  92%  |                                                                              |======================================================================| 100%
rich_pam <- pam2richness(com_pam,which_steps = c(1,5,10))
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
raster::plot(rich_pam)

# }
```
