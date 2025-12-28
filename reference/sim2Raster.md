# sim2Raster: Convert a BAM simulation object to RasterStack

Convert a BAM simulation object to RasterStack.

## Usage

``` r
sim2Raster(sdm_simul, which_steps = NULL, progress_bar = TRUE)
```

## Arguments

- sdm_simul:

  A bam object. See [`sdm_sim`](sdm_sim.md)

- which_steps:

  A numeric vector indicating the simulation steps that are going to be
  converted into raster layers.

- progress_bar:

  Show progress bar

## Value

A RasterStack of species' distribution at each simulation step

## Author

Luis Osorio-Olvera & Jorge Soberón

## Examples

``` r
# \donttest{
model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                          package = "bamm")
model <- raster::raster(model_path)
sparse_mod <- bamm::model2sparse(model,threshold=0.1)
adj_mod <- bamm::adj_mat(sparse_mod,ngbs = 1)
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
#> Simulation progress:
#> [                                                  ][[====>                                             ]   9%[=========>                                        ]  18%[=============>                                    ]  27%[==================>                               ]  36%[======================>                           ]  45%[===========================>                      ]  54%[===============================>                  ]  63%[====================================>             ]  72%[========================================>         ]  81%[=============================================>    ]  90%[==================================================] 100%
sdm_lep_cal_st <- bamm::sim2Raster(sdm_simul = sdm_lep_cal,
                                  which_steps = seq(1,10,by=1))
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%

raster::plot(sdm_lep_cal_st)


# }
```
