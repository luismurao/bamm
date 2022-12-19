## Resubmission

This is a resubmission. In this version I have done what the reviewer asked for:
 
 1) Description: 
     - I changed the Title field to be in title case
     - I changed <doi.org/10.17161/bi.v2i0.4> to <doi:10.17161/bi.v2i0.4>
     - I give more details about what the package does.
     - I added references to the paper describing methods
     - I explained all acronyms in the description text
     - I used single quotes only for software names
 2) Function documentation:
     - All functions have the \value field
     - I gave more details about function return values and added references
 3) Function examples:
     - I replaced dontrun{} with \donttest{}, those examples that take more than 
     5 seconds.
     - I unwrapped examples executable in < 5 sec
 4) Made sure to do not to change the user's options (graphics parameters):
     - I included on.exit functionality to those functions that change 
     default parameter settings for graphics.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release (new submission).

- Added single to some words in the description to avoid a NOTE:
    - Possibly misspelled words in DESCRIPTION: 'Abiotic', 'Biotic' and 'BAM'

## Test environments

- R-hub windows-x86_64-devel (r-devel)
- R-hub  macOS 10.13.6 High Sierra, R-release, CRAN's setup
- windows-latest (release; on GitHub Actions), R 4.2.2
- macOS-latest (release; on GitHub Actions), R 4.2.2
- ubuntu-latest (release; on GitHub Actions), R 4.2.2
- ubuntu-latest (devel; on GitHub Actions), R 4.2.2
- ubuntu-latest (older-1; on GitHub Actions), R 4.2.2
- local ubuntu 20.04.5, install, R 4.1.2
- local OS X 12.4 install, R 4.2.2
- local windows 10 pro install, R 4.2.1

