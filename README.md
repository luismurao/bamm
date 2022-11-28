# The BAM package <a href='https://luismurao.github.io/'><img src='inst/figures/bam_logo.png' align="right" height="139" /></a>

**Authors: Luis Osorio-Olvera & Jorge Soberon**

## Overview

The `bam` package is an R package designed to estimate dynamic models of
species distributions using the concepts of the **BAM** shceme. It allows
to  operate on large matrices (tens of millions of cells) regarding to each 
element of the **BAM**, for example, the adjacency matrix 
(connectivity matrix), and the niche suitability matrices. 

The dynamic model behind the package is the cellular automaton

![\begin{equation} \mathbf{G}_j(t+1) =\mathbf{B}_j(t)\mathbf{A}_j(t) \mathbf{C}_j  \mathbf{G}_j(t) \label{eq:automata} \end{equation}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bequation%7D%20%5Cmathbf%7BG%7D_j(t%2B1)%20%3D%5Cmathbf%7BB%7D_j(t)%5Cmathbf%7BA%7D_j(t)%20%5Cmathbf%7BC%7D_j%20%20%5Cmathbf%7BG%7D_j(t)%20%5Clabel%7Beq%3Aautomata%7D%20%5Cend%7Bequation%7D)

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


```r
if (!require('devtools')) install.packages('devtools')
devtools::install_github('luismurao/bam')
# If you want to build vignette, install pandoc before and then
devtools::install_github('luismurao/bam',build_vignettes=TRUE)
```

## Acknowledgements

We are grateful to our many colleagues in the University of Kansas Niche Modeling
Group for many vivacious and useful discussions on the topics of the paper. 
LOO acknowledges partially supported by Consejo Nacional de Ciencia y Tecnología 
(CONACyT; postdoctoral fellowship number 740751; CVU: 368747).
LOO and JS aknowledges Blitzi Soberon for moral support.
