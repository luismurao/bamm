# permute_pam: Function to permute a Presence-Absence-Matrix.

permute_pam: Function to permute a Presence-Absence-Matrix.

## Usage

``` r
permute_pam(m, niter = NULL, as_sparse = FALSE, randal = "fastball")
```

## Arguments

- m:

  Presence-Absence-Matrix (PAM) or a binary matrix with columns
  representing species and rows sites.

- niter:

  Number of iterations to permute the PAM.

- as_sparse:

  If TRUE the PAM will be returned as a sparse matrix

- randal:

  Randomization algorithm applied to the PAM. Possible choices
  "curveball", "fastball", and "indep_swap".

## Value

Returns a permuted matrix of the same dimensions of m (same number of
rows and columns). Note that the sum of each row and column of this
permuted matrix is equal to that of m. species.

## Details

This function can use the "curveball" (Strona et al., 2014), the
fastball (Godard and Neal, 2022), and the independent swap algorithms.
The implementation of the "fastball" in C++ is provided in
<https://github.com/zpneal/fastball/blob/main/fastball.cpp>. Please when
using the "fastball" algorithm for publications cite Godard and Neal
(2022). When using the "curveball" cite Strona et al. (2014). When using
independent swap ("indep_swap") cite Kembel et al. (2010)

## References

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
ppam <- bamm::permute_pam(m = pam,niter = NULL,as_sparse = FALSE)
# Check if matrices are different
all(pam == ppam)
#> [1] FALSE
# Check if row totals are the same
all(Matrix::rowSums(pam) == Matrix::rowSums(ppam))
#> [1] TRUE
# Check if column total are the same
all(Matrix::colSums(pam) == Matrix::colSums(ppam))
#> [1] TRUE
```
