FROM rocker/tidyverse:4.0.3

# install dependencies
RUN apt-get update && \
  apt-get install -y --no-install-recommends libglpk-dev && \
  apt-get clean

# install CRAN version of dyngen for caching
RUN Rscript -e 'install.packages("dyngen")'

# adding the currently checked out dyngen version
ADD . /home/rstudio/package

# install it
RUN Rscript -e 'install.packages("devtools")' && \
  Rscript -e 'devtools::install_deps("/home/rstudio/package", dependencies = TRUE)' && \
  Rscript -e 'devtools::install("/home/rstudio/package", dependencies = TRUE)' && \
  rm -r /home/rstudio/package

# add getting started script
ADD docker/getting_started.R /home/rstudio/
