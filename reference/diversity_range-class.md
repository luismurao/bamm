# Class `diversity_range`

Class `diversity_range`

## Value

An object of class diversity_range

## Slots

- `alpha`:

  A column vector with species richness per site

- `omega`:

  A column vector with the size of the area of distribution per species.

- `alpha_raster`:

  Species richness in raster format.

- `dispersion_field`:

  A matrix with the set of ranges of all species that occur in at each
  locality.

- `dispersion_field_raster`:

  Raster object with the observed values of dispersion field.

- `diversity_range_raster`:

  Raster object of diversity range.

- `diversity_range_colors`:

  Colors to plot endemism levels.

- `null_dispersion_field_dist`:

  A matrix with dispersion field null distribution.

- `xy_coordinates`:

  A matrix of geographical coordinates

- `n_iterations`:

  Number of iterations used to estimate the dispersion field null
  distribution.

- `nsps`:

  Number of species in the PAM.

- `nsites`:

  Number of sites in the PAM.

## Author

Luis Osorio-Olvera & Jorge Soberón
