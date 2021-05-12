# dyngen 1.0.1

* MINOR CHANGE `.download_cacheable_file()`: Check the return value of `utils::download.file()`, since it is possible that the download will fail with a non-zero status but not an R error. 

* MINOR CHANGE: `kinetics_random_distributions()`: Add function for providing randomised distributions.

## Test environments
* local Fedora install (R release)
* ubuntu 20.04 (with Github Actions; R release, devel)
* Mac OS X (with Github Actions; R release)
* Windows (with Github Actions; R release)
* win-builder (oldrelease, release, devel)

## R CMD check results

── R CMD check results ─────────────────────────────────────── dyngen 1.0.1 ────
Duration: 3m 39.5s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded
