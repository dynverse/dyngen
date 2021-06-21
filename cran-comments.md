# dyngen 1.0.2

* MINOR CHANGE `generate_dataset()`: Fix subplot title.

* NEW FEATURE `plot_summary()`: Create a dedicated function for plotting a dyngen summary.

* BUG FIX `generate_cells()`: Fix incorrect cell count when one of the backbone segments does not have any simulations steps (#26).

* DOCUMENTATION: Update citation to NCOMMS publication.

## Test environments
* local Fedora install (R release)
* ubuntu 20.04 (with Github Actions; R release, devel)
* Mac OS X (with Github Actions; R release)
* Windows (with Github Actions; R release)
* win-builder (oldrelease, release, devel)

## R CMD check results

── R CMD check results ─────────────────────────────────────── dyngen 1.0.2 ────
Duration: 5m 3.1s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded
