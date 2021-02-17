Installation instructions
================

<!-- github markdown built using 
rmarkdown::render("vignettes/installation.Rmd", output_format = rmarkdown::github_document())
-->

## System requirements

This package is supported for Linux, but should also work on Mac OS X
and Windows. It has been tested with [Github
Actions](https://github.com/dynverse/dyngen/actions?query=workflow%3AR-CMD-check)
for R 3.5, 3.6 and 4.0 on the following systems Ubuntu, Windows Server
and Mac OS X.

## Step 1: Installing from CRAN or GitHub (Required)

dyngen is available on CRAN, so you can install it with the following
command.

``` r
install.packages("dyngen")
```

If you would like to install the development version of dyngen from
GitHub instead, run the following command. Use at your own risk!

``` r
install.packages("remotes")
remotes::install_github("dynverse/dyngen@devel", dependencies = TRUE)
```

## Step 2: Configure host system (Recommended)

It’s recommended to let dyngen know where it can cache downloaded files
and how many cores the host system has. If you don’t perform these
steps, running dyngen simulations will take a lot longer than it needs
to.

To do so, start editing your Rprofile by running the following commands:

``` r
install.packages("usethis")
usethis::edit_r_profile()
```

Inside your Rprofile, add the following lines:

``` r
options(Ncpus = 8L) # change this to the number of cores in your system
options(dyngen_download_cache_dir = R_user_dir("dyngen", "data"))
```

After saving and closing this file, restart your R sesstion.

## Step 3: Download cacheable files (Optional)

dyngen will sometimes download files from GitHub and cache them for
later use (if you executed step 2). After installation, you can already
download these files so they will not be downloaded at a later stage.

``` r
library(utils)
library(tools)

zip <- tempfile(fileext = ".zip")
download.file("https://github.com/dynverse/dyngen/archive/data_files.zip", zip)
unzip(zip, exdir = getOption("dyngen_download_cache_dir"), overwrite = TRUE)
```
