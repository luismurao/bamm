# Class for the A (Abiotic) Set of BAM Diagram

Represents the abiotic component (environmental suitability) in the BAM
framework. Contains raster-based niche models and their sparse matrix
representations.

## Value

An object of class setA showClass("setA")

## Slots

- `niche_model`:

  A `RasterLayer` or `RasterStack` representing the environmental
  suitability model (binary or continuous)

- `suit_threshold`:

  Numeric value used to binarize continuous models

- `cellIDs`:

  Numeric vector of raster cell IDs with non-NA values

- `suit_values`:

  Numeric vector of suitability values for continuous models

- `sparse_model`:

  A sparse matrix (`dgCMatrix`) representation of the niche model

## Validity

The `niche_model` slot must be either a `RasterLayer` or `RasterStack`.

## See also

[`model2sparse`](model2sparse.md), [`setM-class`](setM-class.md)

## Author

Luis Osorio-Olvera & Jorge Soberón
