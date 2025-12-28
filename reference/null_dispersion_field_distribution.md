# null_dispersion_field_distribution: Null distribution of the dispersion field

null_dispersion_field_distribution estimates a random distribution of
the dispersion field values.

## Usage

``` r
null_dispersion_field_distribution(
  pam,
  n_iter = 10,
  randal = "fastball",
  parallel = TRUE,
  n_cores = 2
)
```

## Arguments

- pam:

  A Presence-Absence-Matrix of matrix class or sparse matrix.

- n_iter:

  Number of iterations to obtain the distribution.

- randal:

  Randomization algorithm applied to the PAM. Possible choices
  "curveball", "fastball", and "indep_swap".

- parallel:

  If TRUE the computations will be performed in parallel.

- n_cores:

  Number of cores for the parallel computation.

## Value

A data matrix of size nrow(pam) X n_iter with dispersion field values.

## Details

Estimates a random distribution of the dispersion field values. To
obtain random values it uses the function
[`permute_pam`](permute_pam.md) at each step of the iterations.
Randomization of the PAM can be performed using the "fastball" (Godard
and Neal, 2022) and the "curveball" (Strona et al., 2014), and and the
independent swap (Kembel et al. 2010) algorithms. The implementation of
the "fastball" in C++ is provided in
<https://github.com/zpneal/fastball/blob/main/fastball.cpp>

## References

Soberon J, Cavner J (2015). “Indices of Biodiversity Pattern Based on
Presence-Absence Matrices: A GIS Implementation.” *Biodiversity
Informatics*, **10**, 22–34.

Strona G, Nappo D, Boccacci F, Fattorini S, San-Miguel-Ayanz J (2014).
“A fast and unbiased procedure to randomize ecological binary matrices
with fixed row and column totals.” *Nature Communications*, **5**(1),
1–9. ISSN 20411723,
[doi:10.1038/ncomms5114](https://doi.org/10.1038/ncomms5114) ,
<https://www.r-project.org>.

Godard K, Neal ZP (2022). “fastball: a fast algorithm to randomly sample
bipartite graphs with fixed degree sequences.” *Journal of Complex
Networks*, **10**(6), cnac049. ISSN 2051-1329,
[doi:10.1093/comnet/cnac049](https://doi.org/10.1093/comnet/cnac049) ,
https://academic.oup.com/comnet/article-pdf/10/6/cnac049/47758701/cnac049.pdf.

Kembel SW, Cowan PD, Helmus MR, Cornwell WK, Morlon H, Ackerly DD,
Blomberg SP, Webb CO (2010). “Picante: R tools for integrating
phylogenies and ecology.” *Bioinformatics*, **26**, 1463–1464.

## Author

Luis Osorio-Olvera & Jorge Soberón

## Examples

``` r
set.seed(111)
pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
dfield_rand <- bamm::null_dispersion_field_distribution(pam,n_iter=10,
                                                       parallel=FALSE,
                                                       randal="indep_swap",
                                                       n_cores = 2)
#> New names:
#> • `dfield` -> `dfield...1`
#> • `dfield` -> `dfield...2`
#> • `dfield` -> `dfield...3`
#> • `dfield` -> `dfield...4`
#> • `dfield` -> `dfield...5`
#> • `dfield` -> `dfield...6`
#> • `dfield` -> `dfield...7`
#> • `dfield` -> `dfield...8`
#> • `dfield` -> `dfield...9`
#> • `dfield` -> `dfield...10`
head(dfield_rand)
#>      dfrand_1 dfrand_2 dfrand_3 dfrand_4 dfrand_5 dfrand_6 dfrand_7 dfrand_8
#> [1,]      0.6      0.8      0.8      0.8      0.6      0.5      0.5      0.8
#> [2,]      0.8      0.8      1.1      1.0      0.9      1.1      0.8      0.9
#> [3,]      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0
#> [4,]      0.3      0.3      0.3      0.1      0.3      0.3      0.3      0.3
#> [5,]      1.4      1.4      1.2      1.4      1.2      1.4      1.4      1.2
#> [6,]      0.3      0.3      0.3      0.3      0.5      0.5      0.5      0.3
#>      dfrand_9 dfrand_10
#> [1,]      0.6       0.6
#> [2,]      0.7       1.1
#> [3,]      0.0       0.0
#> [4,]      0.5       0.1
#> [5,]      1.4       1.3
#> [6,]      0.5       0.5
```
