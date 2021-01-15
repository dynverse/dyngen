# dyngen 0.4.1

## MINOR CHANGES
* `initialise_model()`: Change defaults of `num_cores` and `download_cache_dir`
  to `getOption("Ncpus")` and `getOption("dyngen_download_cache_dir")` respectively, 
  so you can change the system settings with your R profile.

## BUG FIX
* `wrap_dataset()`: Fix `drop = FALSE` bug when only cell is being sampled.  

## Test environments
* local Fedora 32 install (R 4.0)
* ubuntu 16.04 (with Github Actions; R 3.3, 3.4, 3.5, 3.6, release)
* Mac OS X (with Github Actions; R release)
* Windows (with Github Actions; R release)
* win-builder (oldrelease, release, devel)

## R CMD check results

0 errors | 0 warnings | 0 notes
