# S4 Class Hierarchy for BAM (Biotic-Abiotic-Movement) Modeling

This documentation describes the S4 class system used in the `bamm`
package for modeling species distributions using the BAM
(Biotic-Abiotic-Movement) framework. The classes organize data and
results for ecological niche modeling, dispersal simulation, and
biodiversity analysis.

The foundational class containing geographic coordinates and optional
eigen analysis results for spectral graph methods.

## Slots

- `coordinates`:

  A numeric matrix with two columns (x, y) representing geographic
  coordinates of non-NA cells

- `eigen_vec`:

  Matrix containing eigenvectors from adjacency matrix analysis

- `eigen_val`:

  Numeric vector containing eigenvalues from adjacency matrix analysis

## See also

[`setA-class`](setA.md), [`setM-class`](setM-class.md)

## Author

Luis Osorio-Olvera & Jorge Soberón
