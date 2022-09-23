# dyngen 1.0.5

* MINOR CHANGE: Refactor matrix coercion thanks to Matrix 1.5-0.

* BUG FIX: An update to the dependency package GillespieSSA2 was published to
  attempt to resolve an issue with a malloc free which caused dyngen to be removed
  from CRAN.

## Test environments
* local Fedora install (R release)
* ubuntu 20.04 (with Github Actions; R release, devel)
* Mac OS X (with Github Actions; R release)
* Windows (with Github Actions; R release)
* win-builder (oldrelease, release, devel)

## R CMD check results

── R CMD check results ─────────────────────────────────────── dyngen 1.0.5 ────
Duration: 4m 44.6s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
