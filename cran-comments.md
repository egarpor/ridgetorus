## Test environments

* local R installation, R 4.2.2
* win-builder release
* win-builder devel
* Windows Server 2022, R-release, 32/64 bit
* Windows Server 2022, R-devel, 64 bit
* Windows Server 2022, R-oldrel, 32/64 bit
* Windows Server 2022, R-patched, 32/64 bit
* macOS 10.13.6 High Sierra, R-release, brew
* macOS 10.13.6 High Sierra, R-release, CRAN's setup
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC
* Debian Linux, R-release, GCC
* Debian Linux, R-devel, GCC
* macos-latest (release)
* GH's action R-CMD-check ubuntu-latest (devel)
* GH's action R-CMD-check ubuntu-latest (oldrel-1)
* GH's action R-CMD-check ubuntu-latest (release)
* GH's action R-CMD-check windows-latest (release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.
* NOTE "Possibly misspelled words in DESCRIPTION" has been double checked -- words are correct.
* NOTE "Found the following (possibly) invalid URLs" has been double checked -- URLs are correct.
* Sometimes the examples in ridge_pca() takes slightly more than 10s (like 11s). But this is the main function of the package, so it is useful to have at least one example without \dontrun. The running times on the other functions have been trimmed down to make up for this.
