<!-- badges: start -->
[![R-CMD-check](https://github.com/luismurao/bamm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luismurao/bamm/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/luismurao/bamm/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/luismurao/bamm/actions/workflows/test-coverage.yaml)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/bamm)](https://github.com/r-hub/cranlogs.app)
[![CRAN status](https://www.r-pkg.org/badges/version/bamm)](https://CRAN.R-project.org/package=bamm)
<!-- badges: end -->

# bamm <a href="https://luismurao.github.io/bamm/"><img src="man/figures/logo.png" align="right" height="139" /></a>

**Authors: Luis Osorio-Olvera & Jorge Soberon**

## Overview

The `bamm` package is an R package designed to estimate dynamic models of
species distributions using the concepts of the **BAM** scheme. It allows
to  operate on large matrices (tens of millions of cells) regarding to each 
element of the **BAM**, for example, the adjacency matrix 
(connectivity matrix), and the niche suitability matrices. 

The dynamic model behind the package is the cellular automata

$$\mathbf{G}_j(t+1) =\mathbf{B}_j(t)\mathbf{A}_j(t) \mathbf{C}_j  \mathbf{G}_j(t)$$

The main functions of the package are:

 - *`model2sparse`*: it is the basic function of the package, it converts 
 a binary niche model (in raster format) to sparse matrix model 
 (object of class **setA**). 
    
 - *`adj_mat`*: the function returns the sparse representation of the adjacency 
 matrix of a given raster (generally is the M area but can be any area) given a 
 movement hypothesis. The user can ask the function to return the eigen-analysis 
 of the matrix. 
    
 - *`bam_clusters`*: function to estimate the connectivity of suitable areas 
 given an adjacency matrix. It returns three objects: a) an dynamic map
 (open-street map) of connected areas or clusters; b) the data.frame with 
 coordinates the geographic cluster membership; c) a raster object of with 
 cluster IDs.
     
 - *`csd_estimate`*:This function is used to estimate the CSD-plot. It gives 
 you an idea about the dispersal distance that a species needs to travel to 
 fill its potential area of distribution.
 
 - *`occs2sparse`*: Converts occurrence data into a sparse matrix object. 
 This object is used to declare the initial conditions for modeling the 
 invasion dynamics of a species. 
 
 - *`sdm_sim`*: Simulate single species dispersal dynamics using the cellular 
 automaton of the area of distribution. 
 

# Installation

### CRAN

```r
install.packages("bamm")
```
### GitHub

```r
if (!require('devtools')) install.packages('devtools')
devtools::install_github('luismurao/bamm')
# If you want to build vignette, install pandoc before and then
devtools::install_github('luismurao/bamm',build_vignettes=TRUE)
```

## Acknowledgements 

We are grateful to our many colleagues in the University of Kansas Niche Modeling
Group for many vivacious and useful discussions on the topics of the paper. 
LOO acknowledges partially supported by Programa de Apoyo a Proyectos de 
Investigación e Innovación Tecnológica PAPIIT-IA202824 and Consejo Nacional de 
Ciencia y Tecnología (CONACyT; postdoctoral fellowship number 740751; CVU: 368747).
LOO and JS acknowledges Blitzi Soberon for moral support.
