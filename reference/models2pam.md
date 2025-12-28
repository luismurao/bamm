# models2pam: Converts binary rasters to a PAM

Function to convert binary raster models to a Presence Absences Matrix.

## Usage

``` r
models2pam(
  mods_stack,
  return_coords = FALSE,
  sparse = TRUE,
  parallel = FALSE,
  ncores = 2
)
```

## Arguments

- mods_stack:

  A raster stack containing binary models of each species in the
  community.

- return_coords:

  Logical. If TRUE the pam will be returned with coordinates in the
  first two columns.

- sparse:

  Logical. If TRUE the PAM will be returned as a sparse matrix.

- parallel:

  Logical. If TRUE computations will be done in parallel

- ncores:

  Integer. Number of cores to run the parallel process.

## Value

A presence-absence matrix (PAM).

## Details

For more information about PAM see Soberon and Cavner (2015).

## References

Soberon J, Cavner J (2015). “Indices of Biodiversity Pattern Based on
Presence-Absence Matrices: A GIS Implementation.” *Biodiversity
Informatics*, **10**, 22–34. .

## Author

Luis Osorio-Olvera & Jorge Soberón

## Examples

``` r
# \donttest{
lagos_path <- system.file("extdata/conejos",
                          package = "bamm")
enm_path <- list.files(lagos_path,
                       pattern = ".tif",
                       full.names = TRUE)[1:10]
en_models <- raster::stack(enm_path) >0.01
pam <- bamm::models2pam(en_models,
                        return_coords=TRUE,
                        sparse=FALSE,
                        parallel=FALSE,ncores=1)
head(pam)
#>           x        y Brachylagus_idahoensis_cont Lepus_alleni_cont
#> 1 -91.66667 82.66667                           0                 0
#> 2 -90.66667 82.66667                           0                 0
#> 3 -89.66667 82.66667                           0                 0
#> 4 -88.66667 82.66667                           0                 0
#> 5 -87.66667 82.66667                           0                 0
#> 6 -86.66667 82.66667                           0                 0
#>   Lepus_americanus_cont Lepus_arcticus_cont Lepus_californicus_cont
#> 1                     0                   0                       0
#> 2                     0                   0                       0
#> 3                     0                   0                       0
#> 4                     0                   0                       0
#> 5                     0                   0                       0
#> 6                     0                   0                       0
#>   Lepus_callotis_cont Lepus_europaeus_cont Lepus_othus_cont
#> 1                   0                    0                0
#> 2                   0                    0                0
#> 3                   0                    0                0
#> 4                   0                    0                0
#> 5                   0                    0                0
#> 6                   0                    0                0
#>   Lepus_townsendii_cont Ochotona_collaris_cont
#> 1                     0                      0
#> 2                     0                      0
#> 3                     0                      0
#> 4                     0                      0
#> 5                     0                      0
#> 6                     0                      0
# }
```
