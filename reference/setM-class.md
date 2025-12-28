# Class for the M set of the `bamm` diagram

Class for the M set of the `bamm` diagram

## Value

An object of class setM

## Slots

- `adj_matix`:

  An adjacency matrix

- `adj_list`:

  An adjacency list

- `initial_points`:

  A presence-absence vector with species' occurrences

- `n_initial_points`:

  Number of initial points used to start the dispersal process

- `ngbs`:

  Number of neighbors

## Author

Luis Osorio-Olvera & Jorge Soberón

## Examples

``` r
showClass("setM")
#> Class "setM" [package "bamm"]
#> 
#> Slots:
#>                                                                   
#> Name:      adj_matrix       adj_list initial_points           ngbs
#> Class:      dgCMatrix           list      dgCMatrix        numeric
#>                                                    
#> Name:     coordinates      eigen_vec      eigen_val
#> Class:         matrix         matrix        numeric
#> 
#> Extends: "g_area"
#> 
#> Known Subclasses: "bam"
```
