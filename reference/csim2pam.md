# csim2pam: Converts community simulation to a Presence Absence Matrix (PAM)

Converts community simulation object into a Presence Absence Matrices
(PAM) for a given simulation steps.

## Usage

``` r
csim2pam(community_sim, which_steps, progress_bar = TRUE)
```

## Arguments

- community_sim:

  An object of class [`community_bam`](community_sim-class.md).

- which_steps:

  Steps in the simulation object to be converted into a PAM

- progress_bar:

  Show progress bar

## Value

An object of class [`pam`](pam-class.md); it contains five slots. 1)
pams: a list of sparse matrices with Presence-Absence information
(PAMs). 2) which_steps: time steps corresponding to each PAM. 3)
sp_names: a vector of species names. 4) the grid area used in the
simulation. 5) Non NA cell (pixel) IDs.

## Details

For details about the object community_sim see
[`community_sim`](community_sim.md)

## References

Soberón J, Osorio-Olvera L (2023). “A dynamic theory of the area of
distribution.” *Journal of Biogeography6*, **50**, 1037-1048.
[doi:10.1111/jbi.14587](https://doi.org/10.1111/jbi.14587) ,
<https://onlinelibrary.wiley.com/doi/abs/10.1111/jbi.14587>. .

## Author

Luis Osorio-Olvera & Jorge Soberón

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
nsteps <- 10
sdm_comm <- bamm::community_sim(en_models = en_models,
                               ngbs_vect = ngbs_vect,
                               init_coords = init_coords,
                               nsteps = nsteps,
                               threshold = 0.1)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%

pamt10 <- bamm::csim2pam(community_sim = sdm_comm ,
                        which_steps = 10)
#>   |                                                                              |                                                                      |   0%
pams <- bamm::csim2pam(community_sim = sdm_comm ,
                       which_steps = seq_len(10))
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
rich_pam <- bamm::pam2richness(pams,which_steps = c(1,5))
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
print(rich_pam)
#> class      : RasterStack 
#> dimensions : 75, 112, 8400, 2  (nrow, ncol, ncell, nlayers)
#> resolution : 1, 1  (x, y)
#> extent     : -168.1667, -56.16667, 8.166667, 83.16667  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +no_defs 
#> names      : time_step_1, time_step_5 
#> min values :           0,           0 
#> max values :           1,           3 
#> 
# }
```
