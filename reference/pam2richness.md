# pam2richness: Converts Presence Absence Matrix (pam object) to richness raster

Converts Presence Absence Matrix (pam object) to richness raster

## Usage

``` r
pam2richness(pamobj, which_steps)
```

## Arguments

- pamobj:

  An object of class pam see [`csim2pam`](csim2pam.md)

- which_steps:

  Time steps in the pam to convert

## Value

A RasterStack richness for each simulation step

## Author

Luis Osorio-Olvera & Jorge Soberón.

## Examples

``` r
# \donttest{
lagos_path <- system.file("extdata/conejos",
                          package = "bamm")
enm_path <- list.files(lagos_path,
                       pattern = ".tif",
                    full.names = TRUE)[seq(1,10)]
en_models <- raster::stack(enm_path)
ngbs_vect <- sample(2,replace = TRUE,
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

pams <-bamm::csim2pam(community_sim = sdm_comm ,
                      which_steps = seq_len(nsteps))
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
richness_stack <- bamm::pam2richness(pamobj = pams,
                                     which_steps = pams@which_steps)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
raster::plot(richness_stack)

# }
```
