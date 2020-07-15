FROM rocker/tidyverse:4.0.2

ARG GITHUB_PAT

RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")' && \
  # install dyngen
  Rscript -e 'remotes::install_github("dynverse/dyngen@devel", dependencies = TRUE, repos = "https://cran.rstudio.com")' && \
  # install dyno and slingshot (for example.R)
  Rscript -e 'remotes::install_github(c("dynverse/ti_slingshot/package", "dynverse/dyno@devel"), repos = "https://cran.rstudio.com")' && \
  # set default dynwrap backend to r_wrapper (for example.R)
  echo "options(dynwrap_backend = 'r_wrapper')" >> /usr/local/lib/R/etc/Rprofile.site
  
ADD example.R /home/rstudio/
