# dyngen 1.1.0

* BREAKING CHANGE: Minimum R version bumped to 4.1.0.

* MINOR CHANGE `as_anndata()`: Replace the `anndata` R package (Python-based)
  with `anndataR` (Bioconductor), removing the Python dependency from the
  `anndata` output format.

* MINOR CHANGE `plot_simulation_expression()`: Move `ggrepel` from `Imports`
  to `Suggests`.

* MINOR CHANGE: Replace `magrittr::%>%` with the native pipe `|>`, deprecated
  dplyr/tidyr verbs with their modern equivalents, and deprecated ggplot2 APIs
  with current ones throughout the package.

* BUG FIX `get_timings()`: Compatibility with dplyr >= 1.1.0.

* BUG FIX `bblego_linear()`: Fix broken pipe chain in `doublerep2` type.

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
