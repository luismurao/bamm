# jaccard: Estimates the Jaccard index for comparing two binary maps

Estimates the Jaccard index for comparing two binary maps

## Usage

``` r
jaccard(m1, m2)
```

## Arguments

- m1:

  A binary raster A or an object of class setA returned by the function
  [`model2sparse`](model2sparse.md).

- m2:

  A binary raster A or an object of class setA returned by the function
  [`model2sparse`](model2sparse.md).

## Value

Returns a data.frame with three values: 1) jaccard (Jaccard index), 2)
percentage_m1 (the percentage of m1 that the intersection \\\|A \cap
B\|\\ represents), and 3) percentage_m2

## Details

The Jaccard index is computed as follows \$\$J(A,B)={{\|A\cap
B\|}\over{\|A\cup B\|}}={{\|A\cap B\|}\over{\|A\|+\|B\|-\|A\cap B\|}}.
\$\$

## Author

Luis Osorio-Olvera & Jorge Soberón

## Examples

``` r
m1_path <- system.file("extdata/conejos/Lepus_othus_cont.tif",
                       package = "bamm")
m2_path <- system.file("extdata/conejos/Brachylagus_idahoensis_cont.tif",
                       package = "bamm")
m1 <- raster::raster(m1_path) > 0.01
m2 <- raster::raster(m2_path) >0.01
jcc <- bamm::jaccard(m1,m2)
print(jcc)
#>      jaccard percentage_of_m1 percentage_of_m2
#> 1 0.04545455         5.475763         21.10727
```
