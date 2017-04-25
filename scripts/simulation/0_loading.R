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

library(dambiutils)

library(dyngen)


devtools::load_all(".")