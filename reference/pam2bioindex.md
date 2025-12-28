# pam2bioindex: PAM to biodiversity index

pam2bioindex estimates various biodiversity indices for a certain PAM.

## Usage

``` r
pam2bioindex(pam, biodiv_index = "dispersion_field", as_sparse = FALSE)
```

## Arguments

- pam:

  A Presence-Absence-Matrix of matrix class or sparse matrix.

- biodiv_index:

  Possible values are alpha, omega, wbeta (Whittaker’s multiplicative
  beta index), laBeta (Lande’s additive beta index) dispersion_field,
  all.

- as_sparse:

  Return indices as sparse objects

## Value

An object of class [`bioindex`](bioindex-class.md) with three slots each
represents a matrix of diversity indices: alpha, omega, and dispersion
field, richness_field.

## Details

The biodiversity indices can be found in Soberón and Cavner (2015).

## References

Soberon J, Cavner J (2015). “Indices of Biodiversity Pattern Based on
Presence-Absence Matrices: A GIS Implementation.” *Biodiversity
Informatics*, **10**, 22–34.

## Author

Luis Osorio-Olvera & Jorge Soberón

## Examples

``` r
set.seed(111)
pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
bioindices <- bamm::pam2bioindex(pam=pam,biodiv_index="all")
# Return results as sparse models
bioindices <- bamm::pam2bioindex(pam=pam,biodiv_index="all",as_sparse=TRUE)
bioindices@alpha
#> 10 x 1 Matrix of class "dgeMatrix"
#>       [,1]
#>  [1,]    2
#>  [2,]    3
#>  [3,]    0
#>  [4,]    1
#>  [5,]    4
#>  [6,]    1
#>  [7,]    1
#>  [8,]    5
#>  [9,]    3
#> [10,]    5
bioindices@omega
#> 10 x 1 Matrix of class "dgeMatrix"
#>       [,1]
#>  [1,]    1
#>  [2,]    1
#>  [3,]    2
#>  [4,]    3
#>  [5,]    3
#>  [6,]    3
#>  [7,]    1
#>  [8,]    3
#>  [9,]    3
#> [10,]    5
bioindices@dispersion_field
#> 10 x 1 Matrix of class "dgeMatrix"
#>       [,1]
#>  [1,]  0.6
#>  [2,]  0.9
#>  [3,]  0.0
#>  [4,]  0.3
#>  [5,]  1.1
#>  [6,]  0.3
#>  [7,]  0.5
#>  [8,]  1.3
#>  [9,]  1.0
#> [10,]  1.7
```
