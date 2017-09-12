source("https://bioconductor.org/biocLite.R")
library(scater) # stupid package overrides filter, mutate and other useful function, WHY???? Nobody knows!
library(tidyverse)
library(dplyr)
library(fastgssa)
library(stringr)
library(pheatmap)
library(parallel)
library(PRISM)
# qsub_config <- PRISM::create_qsub_config(
#   remote = "prism",
#   local_tmp_path = "/home/wouters/thesis/tmp/.r2gridengine",
#   remote_tmp_path = "/scratch/irc/personal/wouters/.r2gridengine",
#   execute_before = c("module unload python", "module load python"), 
#   memory = "4G"
# )
qsub_config <- PRISM::override_qsub_config(
  memory = "4G"
)
PRISM::set_default_qsub_config(qsub_config, permanent = F)

library(readr)
library(Biobase)
library(igraph)
library(cowplot)
library(magrittr)
library(purrr)
library(pdist)

library(dambiutils)


devtools::load_all(".")





# install on prism devtools::install_github("Zouter/dyngen", auth_token="1381e4efbbcc986c4601f8a457943f1b16e31f58")
# devtools::install_github("Zouter/dyneval", auth_token="1381e4efbbcc986c4601f8a457943f1b16e31f58")

# devtools::install_github("dambi/fastgssa", host="github.ugent.be/api/v3", auth_token="8aaf33b05a1d0c728f070897ca344f468d17fcfe")
# devtools::install_github("dambi/dambiutils", host="github.ugent.be/api/v3", auth_token="8aaf33b05a1d0c728f070897ca344f468d17fcfe")

# devtools::install_github("rcannood/PRISM", auth_token="1381e4efbbcc986c4601f8a457943f1b16e31f58")