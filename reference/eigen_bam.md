# eigen_bam: Compute the Eigen system of two bam objects

Calculates the Eigen values and Eigen vectors of bam objects

## Usage

``` r
eigen_bam(A = NULL, M = NULL, AM = TRUE, which_eigen = 1, rmap = TRUE)
```

## Arguments

- A:

  A bam object of class setA.

- M:

  A bam object of class setM.

- AM:

  A logical value to specify whether to use the product AM or MA. If
  true the AM will be returned else the product MA will be returned.

- which_eigen:

  An integer representing the which eigen value and eigen vector will be
  computed.

- rmap:

  Logical. If TRUE the function will return a map of the eigen vector of
  the product AM.

## Value

A list with four objects. 1) eigen_values (these are indicated in
which_eigen parameter of the function), 2) eigen_vectors (the
corresponding eigen vectors of each eigen value), 3) Standardized eigen
vectors (0 to 1), 4) A RasterLayer depicting the information of the
first eigen vector of the system.

## Details

The eigenvector associated with the dominant eigenvalue of an adjacency
matrix provides information about the number of forms in which a cell
can be visited from other cells. Details about the eigen analysis in the
context of the area of distribution can be found in Soberon and
Osorio-Olvera (2022).

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
model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                          package = "bamm")
model <- raster::raster(model_path)
sparse_mod <- bamm::model2sparse(model = model,0.75)
plot(sparse_mod@niche_model)

adj_mod <- bamm::adj_mat(sparse_mod,ngbs = 1,eigen_sys = TRUE)
# Product AM
eig_bam_am <- bamm::eigen_bam(A=sparse_mod,M=adj_mod,AM=TRUE)
raster::plot(eig_bam_am$map)

# Product MA
eig_bam_ma <- bamm::eigen_bam(A=sparse_mod,M=adj_mod,AM=FALSE)
raster::plot(eig_bam_ma$map)

# }
```
