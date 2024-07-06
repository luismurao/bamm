# bamm 0.5.0

  - The new version of the package implements the fastball and the independent 
    swap algorithms to randomize the presence-absence matrix while preserving 
    fixed richness and species incidence. The fastball and the independent 
    swap algorithms run 10 and 60 times faster than the package implementation 
    of the curveball algorithm. This is particularly relevant for running the
    the null_dispersion_field and the diversity_range_analysis functions.
    
 - The pam2bioindex function now allows the estimation of nestedness index

# bamm 0.4.2.6

* Added a `NEWS.md` file to track changes to the package.

* Tested in:
    - R-hub windows-x86_64-devel (r-devel)
    - R-hub  macOS 10.13.6 High Sierra, R-release, CRAN setup
    - windows-latest (release; on GitHub Actions), R 4.2.2
    - macOS-latest (release; on GitHub Actions), R 4.2.2
    - ubuntu-latest (release; on GitHub Actions), R 4.2.2
    - ubuntu-latest (devel; on GitHub Actions), R 4.2.2
    - ubuntu-latest (older-1; on GitHub Actions), R 4.2.2
    - local ubuntu 20.04.5, install, R 4.1.2
    - local OS X 12.4 install, R 4.2.2
* Package website build using pkgdown:
    - https://luismurao.github.io/bamm
