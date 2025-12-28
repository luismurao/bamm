# occs2sparse: Converts occurrence data into a sparse matrix object

occs2sparse: Converts occurrence data into a sparse matrix object

## Usage

``` r
occs2sparse(modelsparse, occs)
```

## Arguments

- modelsparse:

  A setA object returned by the function
  [`model2sparse`](model2sparse.md)

- occs:

  A matrix or a data.frame containing two columns. The first one is the
  longitude and the second is the latitude.

## Value

A sparse vector of zeros (presences) and ones (absences).

## Details

Rows of this column vector represent non NA pixels of the niche model.

## Author

Luis Osorio-Olvera & Jorge Soberón

## Examples

``` r
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

head(occs_sparse)
#> 6 x 1 sparse Matrix of class "dgCMatrix"
#>       
#> [1,] .
#> [2,] .
#> [3,] .
#> [4,] .
#> [5,] .
#> [6,] .
```
