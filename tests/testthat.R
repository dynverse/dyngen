library(testthat)
library(dyngen)
library(purrr)

Sys.setenv("R_TESTS" = "")

test_check("dyngen")

