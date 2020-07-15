FROM rocker/tidyverse:4.0.2

ARG GITHUB_PAT

RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")' && \
  Rscript -e 'remotes::install_github("dynverse/dyngen@devel", dependencies = TRUE, repos = "https://cran.rstudio.com")'
  
RUN  Rscript -e 'remotes::install_github("dynverse/dyngen@devel", dependencies = TRUE, repos = "https://cran.rstudio.com")'
  
ADD example.R /home/rstudio/
