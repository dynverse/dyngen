library(tidyverse)
library(PRISM)
library(dynalysis)
#library(dyngen)

dataset_id <- "synthetic/v7"
dataset_preprocessing(dataset_id)
remote_folder <- "/group/irc/shared/dynalysis/" # dynalysis on prism

unlink(dataset_file(), recursive=TRUE);dir.create(dataset_file(), recursive=TRUE, showWarnings = FALSE)
unlink(dataset_preproc_file(), recursive=TRUE);dir.create(dataset_preproc_file(), recursive=TRUE, showWarnings = FALSE)
PRISM:::run_remote(glue::glue("rm -r {paste0(remote_folder, dataset_file(relative = TRUE))}"), "prism")
PRISM:::run_remote(glue::glue("rm -r {paste0(remote_folder, dataset_preproc_file(relative = TRUE))}"), "prism")
PRISM:::run_remote(glue::glue("mkdir {paste0(remote_folder, dataset_preproc_file(relative = TRUE))}"), "prism")
PRISM:::run_remote(glue::glue("mkdir {paste0(remote_folder, dataset_file(relative = TRUE))}"), "prism")

# Prepare environment for remote -----------------
ncores <- 6
prepare_environment <- function(load = character(), params_i = 1) {
  library(tidyverse)
  library(dynalysis)
  
  options(ncores=6)
  
  dataset_preprocessing(dataset_id)
}
prepare_environment()

load_data <- function(load = character(), params_i = 1) {
  map(load, function(id) {
    readRDS(dataset_preproc_file(pritt("{params_i}_{id}.rds")))
  }) %>% set_names(load)
}

# Generate param settings ----------------------------
settings_model <- tribble(
  ~modulenet_name, ~totaltime,
  "linear", 5,
  "bifurcating", 5,
  "linear_long", 30,
  "cycle", 30,
  "consecutive_bifurcating", 10,
  "bifurcating_converging", 15,
  "trifurcating", 10,
  "converging", 10,
  "bifurcating_loop", 30
)
settings_platform <- tibble(platform_name = list.files("ext_data/platforms/") %>% gsub("(.*)\\.rds", "\\1", .)) %>% filter(row_number() <= 10) %>% bind_rows(tibble(platform_name = "small"), .)

settings_replicates <- tibble(replicate_id = 1)
settings <- tidyr::crossing(settings_model, settings_platform, settings_replicates)
settings <- settings %>% 
  group_by(modulenet_name) %>% 
  mutate(dataset_id = paste0(modulenet_name, "_", seq_len(n()))) %>% 
  ungroup() %>% 
  mutate(params_i = seq_len(n()))

settings <- settings[1, ]

update_params <- function(base_params=dyngen:::base_params, ...) {
  dots <- list(...)
  
  if("modulenet_name" %in% names(dots)) base_params$model$modulenet_name <- dots$modulenet_name
  if("totaltime" %in% names(dots)) base_params$simulation$totaltime <- dots$totaltime
  if("platform_name" %in% names(dots)) base_params$model$platform <- readRDS(paste0(find.package("dyngen"), "/ext_data/platforms/", dots$platform_name, ".rds"))
  if("platform_name" %in% names(dots)) base_params$experiment$platform <- readRDS(paste0(find.package("dyngen"), "/ext_data/platforms/", dots$platform_name, ".rds"))
  
  base_params$settings <- dots
  
  base_params
}

paramsets <- map(seq_len(nrow(settings)), function(row_id) {
  row <- dynutils::extract_row_to_list(settings, row_id)
  invoke(update_params, row)
})

saveRDS(paramsets, dataset_preproc_file("paramsets.rds"))
paramsets <- readRDS(dataset_preproc_file("paramsets.rds"))
PRISM:::rsync_remote("", dataset_preproc_file(), "prism", paste0(remote_folder, dataset_preproc_file(relative=TRUE)))

params <- paramsets[[1]]
params_i = params$settings$params_i






## Run on cluster ------
qsub_environment <- list2env(lst(ncores, prepare_environment, dataset_id))

qsub_config <- override_qsub_config(
  num_cores = ncores, 
  memory = paste0("10G"), 
  wait=FALSE, 
  r_module=NULL, 
  execute_before="", 
  name = "dyngen", 
  stop_on_error = F, 
  max_wall_time = "24:00:00"
)
qsub_config_single <- override_qsub_config(qsub_config, num_cores = 1)
qsub_packages <- c("tidyverse", "dynalysis", "dyngen")

run_cluster <- function(func, output_file, qsub_config=qsub_config, paramsets) {
  wrapper <- function(params) {
    if (!file.exists(dataset_preproc_file(pritt("{params$settings$params_i}_{output_file}")))) {
      func(params)
    } else {
      TRUE
    }
  }
  
  output_name <- gsub("(.*)\\..*", "\\1", output_file)
  
  handle <- qsub_lapply(
    qsub_config=qsub_config %>% list_modify(name = output_name), 
    qsub_environment=qsub_environment, 
    qsub_packages=qsub_packages,
    paramsets,
    func
  )
  handle
}


model_handle <- run_cluster(generate_model, "model.rds", qsub_config_single, paramsets)
simulate_handle <- run_cluster(simulate, "simulate.rds", qsub_config, paramsets)
check_simulation_handle <- run_cluster(check_simulation, "simulation_checks.rds", qsub_config_single %>% list_modify(memory="15G"), paramsets)
gs_handle <- run_cluster(extract_goldstandard, "gs.rds", qsub_config_single %>% list_modify(memory="30G"), paramsets)
run_cluster(generate_experiment, "experiment.rds", qsub_config_single %>% list_modify(memory="20G"), paramsets)
run_cluster(normalise, "normalisation.rds", qsub_config, paramsets)
run_cluster(plot_goldstandard, "gs_plot.rds", qsub_config_single %>% list_modify(memory="20G"), paramsets)
run_cluster(plot_experiment, "experiment_plot.rds", qsub_config_single, paramsets)
run_cluster(wrap, "task.rds", qsub_config_single %>% list_modify(memory="20G"), paramsets)

qsub_retrieve(simulate_handle)




## Group tasks
tasks <- list.files(dataset_preproc_file(), "*task.rds", full.names=T) %>% sort() %>% map(readRDS) %>% dynutils::list_as_tibble()
tasks <- tasks %>% filter(map(settings, "platform_name") != "small")

write_rds(tasks, dataset_file("tasks.rds"))













handle <- qsub_lapply(
  qsub_config=qsub_config_single %>% list_modify(name = "model"), 
  qsub_environment=qsub_environment, 
  qsub_packages=qsub_packages, 
  seq_along(paramsets),
  function(params_i) {
    if (!file.exists(pritt("{params_id}_model.rds"))) {
      params <- paramsets[[params_i]]
      options(ncores = 1)
    
      model <- invoke(dyngen:::generate_model_from_modulenet, params$model)
      model$uuids <- list(model=uuid::UUIDgenerate())
      saveRDS(model, dataset_preproc_file(pritt("{params_id}_model.rds")))
    }
  }
) %>% saveRDS("model_handle.rds")
models <- qsub_retrieve(readRDS("model_handle.rds"))
PRISM:::rsync_remote("prism", remote_folder, "", folder)

## SIMULATE CELLS ---------------------------
# walk(seq_along(paramsets), function(params_i) {
qsub_lapply(
  qsub_config = qsub_config %>% list_modify(name = "simulation"), 
  qsub_environment=qsub_environment, 
  qsub_packages = qsub_packages, 
  seq_along(paramsets),
  function(params_i) {
    if (!file.exists(pritt("{params_id}_simulation.rds"))) {
      params <- paramsets[[params_i]]
      model <- readRDS(pritt("{params_id}_model.rds"))
      options(ncores = ncores)
  
      simulation <- invoke(dyngen:::simulate_multiple, params$simulation, model$system)
      simulation$uuids <- list(simulation=uuid::UUIDgenerate()) %>% c(model$uuids)
      saveRDS(simulation, pritt("{params_id}_simulation.rds"))
    }
    TRUE
  }
) %>% saveRDS("simulations_handle.rds")
simulations <- qsub_retrieve(readRDS("simulations_handle.rds"))
PRISM:::rsync_remote("prism", remote_folder, "", folder)

## Check simulations -------
PRISM::qsub_run(qsub_config = qsub_config %>% list_modify(name = "simulation_check", memory="15G"), qsub_environment=qsub_environment, qsub_packages = qsub_packages, function(i) {
  parallel::mclapply(seq_along(paramsets), function(params_i) {
    params <- paramsets[[params_i]]
    data <- list(exists = FALSE)
    if (file.exists(pritt("{params_id}_simulation.rds"))) {
      data$exists <- TRUE
      simulation <- readRDS(pritt("{params_id}_simulation.rds"))
      tibble(
        nsteps = nrow(simulation$stepinfo),
        params_i = params_i
      )
    }
  }, mc.cores=ncores)
}) %>% saveRDS("simulations_checks_handle.rds")
checks <- qsub_retrieve(readRDS("simulations_checks_handle.rds"))


# GOLD STANDARD ----------------------
# walk(seq_along(paramsets), function(params_i) {
qsub_lapply(
  qsub_config = qsub_config_single %>% list_modify(name = "gs", memory="30G"), 
  qsub_environment=qsub_environment, 
  qsub_packages = qsub_packages, 
  seq_along(paramsets), 
  function(params_i) {
    if (!file.exists(pritt("{params_id}_gs.rds"))) {
      params <- paramsets[[params_i]]
      model <- readRDS(pritt("{params_id}_model.rds"))
      simulation <- readRDS(pritt("{params_id}_simulation.rds"))
      options(ncores = 1)
    
      print("Preprocessing")
      simulation <- dyngen:::preprocess_simulation_for_gs(simulation, model, params$gs$smooth_window) # do preprocessing separate, otherwise zoo will stay in an infinite loop in case of later error
      gs <- invoke(dyngen:::extract_goldstandard, params$gs, simulation, model, preprocess=FALSE)
      gs$checks <- dyngen:::check_goldstandard(gs)
      gs$uuids <- list(gs=uuid::UUIDgenerate()) %>% c(simulation$uuids)
      saveRDS(gs, pritt("{params_id}_gs.rds"))
    }
  }
) %>% saveRDS("gs_handle.rds")
gs <- qsub_retrieve(readRDS("gs_handle.rds"))
PRISM:::rsync_remote("prism", remote_folder, "", folder)

checks <- map_dfr(seq_along(paramsets), function(params_i) {
  print(glue::glue("{params_i} / {length(paramsets)} ======================================"))
  gs <- readRDS(pritt("{params_id}_gs.rds"))
  
  tibble(all_represented = all(gs$checks$edge_counts > 0))
})

# EXPERIMENT ----------------------------------
# walk(seq_along(paramsets), function(params_i) {
handle <- qsub_lapply(qsub_config = qsub_config_single %>% list_modify(name = "experiment", memory="20G"), qsub_environment=qsub_environment, seq_along(paramsets), function(params_i) {
  if (!file.exists(pritt("{params_id}_experiment.rds"))) {
    library(tidyverse);library(dyngen)
    params <- paramsets[[params_i]]
    params$experiment %>% list2env(.GlobalEnv)
    
    simulation <- readRDS(pritt("{params_id}_simulation.rds"))
    gs <- readRDS(pritt("{params_id}_gs.rds"))
    options(ncores = 1)
    
    experiment <- invoke(dyngen:::run_experiment, params$experiment, simulation, gs)
    experiment$uuids <- list(experiment=uuid::UUIDgenerate()) %>% c(gs$uuids)
    saveRDS(experiment, pritt("{params_id}_experiment.rds"))
  }
}) %>% saveRDS("experiments_handle.rds")
experiments <- qsub_retrieve(readRDS("experiments_handle.rds"))

# NORMALISE ----------------------------
qsub_lapply(qsub_config = qsub_config_single %>% list_modify(name = "normalisat"), qsub_environment=qsub_environment, qsub_packages = qsub_packages, seq_along(paramsets), function(params_i) {
  if (!file.exists(pritt("{params_id}_normalisation.rds"))) {
    params <- paramsets[[params_i]]
    experiment <- readRDS(pritt("{params_id}_experiment.rds"))
    options(ncores = 1)
    
    graphics.off()
    pdf(pritt("{params_id}_normalisation_plot.pdf"), width=12, height=12)
    dev.control('enable')
    normalisation <- invoke(dynutils::normalise_filter_counts, params$normalisation, experiment$counts, verbose = TRUE)
    normalisation$uuids <- list(normalisation=uuid::UUIDgenerate()) %>% c(experiment$uuids)
    saveRDS(normalisation, pritt("{params_id}_normalisation.rds"))
    
    walk(normalisation$normalisation_plots, print)
    dyngen:::plot_normalisation(experiment, normalisation)
    graphics.off()
  }
}) %>% saveRDS("normalisation_handle.rds")
normalizations <- qsub_retrieve(readRDS("normalisation_handle.rds"))

# Plot GOLD STANDARD ----------------------------------------
qsub_lapply(qsub_config = qsub_config_single %>% list_modify(name = "gs_plot", memory="15G"), qsub_environment=qsub_environment, qsub_packages = qsub_packages, seq_along(paramsets), function(params_i) {
  tryCatch({
    if (!file.exists(pritt("{params_id}_gs_plot.pdf"))) {
      params <- paramsets[[params_i]]
      model <- readRDS(pritt("{params_id}_model.rds"))
      simulation <- readRDS(pritt("{params_id}_simulation.rds"))
      gs <- readRDS(pritt("{params_id}_gs.rds"))
      
      graphics.off()
      pdf(pritt("{params_id}_gs_plot.pdf"), width=12, height=12)
      dyngen:::plot_net(model, label=FALSE, main_only = FALSE)
      dyngen:::plot_modulenet(model)
      dyngen:::plot_goldstandard(simulation, model, gs)
      graphics.off()
    }
  }, error=function(e) {print(params_i)}, finally={graphics.off()})
}) %>% saveRDS("gs_plot_handle.rds")
gs <- qsub_retrieve(readRDS("gs_plot_handle.rds"))

# Plot EXPERIMENT ----------------------------------------
params_i <- 3
qsub_lapply(qsub_config = qsub_config_single %>% list_modify(name = "exp_plot"), qsub_environment=qsub_environment, qsub_packages = qsub_packages, seq_along(paramsets), function(params_i) {
  tryCatch({
    params <- paramsets[[params_i]]
    experiment <- readRDS(pritt("{params_id}_experiment.rds"))
    
    graphics.off()
    pdf(pritt("{params_id}_experiment_plot.pdf"), width=12, height=12)
    dyngen:::plot_experiment(experiment)
    graphics.off()
  }, error=function(e) {print(params_i)}, finally={graphics.off()})
}) %>% saveRDS("experiment_plot_handle.rds")
experiment_plot <- qsub_retrieve(readRDS("experiment_plot_handle.rds"))

# Wrap into tasks ----------------------------
qsub_lapply(qsub_config = qsub_config_single %>% list_modify(name = "wrap", memory="20G"), qsub_environment=qsub_environment, qsub_packages = qsub_packages, seq_along(paramsets), function(params_i) {
  if (!file.exists(pritt("{params_id}_task.rds"))) {
    params <- paramsets[[params_i]]
    model <- readRDS(pritt("{params_id}_model.rds"))
    simulation <- readRDS(pritt("{params_id}_simulation.rds"))
    gs <- readRDS(pritt("{params_id}_gs.rds"))
    experiment <- readRDS(pritt("{params_id}_experiment.rds"))
    normalisation <- readRDS(pritt("{params_id}_normalisation.rds"))
    
    task <- dyngen::wrap_task(params, model, simulation, gs, experiment, normalisation)
    task$uuids <- list(task=uuid::UUIDgenerate()) %>% c(normalisation$uuids)
    
    task %>% saveRDS(pritt("{params_id}_task.rds"))
  }
}) %>% saveRDS("wrap_handle.rds")
wrap <- qsub_retrieve(readRDS("wrap_handle.rds"))

tasks <- list.files(folder, "*task.rds", full.names=T) %>% sort() %>% map(readRDS) %>% dynutils::list_as_tibble()
tasks <- tasks %>% filter(map(info, "platform_name") != "small")

write_rds(tasks, dataset_file("tasks.rds"))
tasks <- read_rds("../dynalysis/analysis/data/derived_data/datasets/synthetic/v6.rds")

tasks$trajectory_type %>% table()
