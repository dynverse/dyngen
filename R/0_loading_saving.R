#' @export
save_model = function(model, modelid=model$info$id) {
  model_folder = file.path(.datasets_location, "models/", modelid)
  dir.create(model_folder)
  saveRDS(model, file.path(model_folder, "model.rds"))
  
  models = readRDS(file.path(.datasets_location, "models.rds"))
  models %>% filter(id != modelid) %>% bind_rows(model$info) %>% saveRDS(file.path(.datasets_location, "models.rds"))
}
#' @export
load_model = function(modelid, contents=contents_model()) {
  model_folder = file.path(.datasets_location, "models/", modelid)
  readRDS(file.path(model_folder, "model.rds"))
}

#' @export
contents_experiment = function(molecules=FALSE, cellinfo=TRUE, expression=TRUE, simulations=FALSE, model=TRUE, expression_modules=FALSE, info=TRUE) {
  as.list(environment(), all=TRUE)
}
#' @export
load_experiment = function(experimentid, contents=contents_experiment()) {
  experiment = list()
  experiment_folder = file.path(.datasets_location, "/experiments/", experimentid)
  loadin = names(contents)[purrr::map_lgl(contents, ~.)]
  experiment = purrr::map(loadin, ~readRDS(file.path(experiment_folder, paste0(., ".rds"))))
  names(experiment) = loadin
  experiment
}
#' @export
save_experiment = function(experiment, experimentid=experiment$info$id) {
  experiment_folder = file.path(.datasets_location, "experiments/", experimentid)
  dir.create(experiment_folder)
  purrr::walk(names(experiment), ~saveRDS(experiment[[.]], file.path(experiment_folder, paste0(., ".rds"))))
  
  experiments = readRDS(file.path(.datasets_location, "experiments.rds"))
  experiments %>% filter(id != experimentid) %>% bind_rows(experiment$info) %>% saveRDS(file.path(.datasets_location, "experiments.rds"))
  print(experiments %>% filter(id != experimentid) %>% bind_rows(experiment$info))
}
#' @export
delete_experiments = function(experimentids) {
  purrr::walk(experimentids, ~unlink(file.path(.datasets_location, "experiments/", .), T))
  experiments = readRDS(file.path(.datasets_location, "experiments.rds"))
  experiments %>% filter(!(id %in% experimentids))  %>% saveRDS(file.path(.datasets_location, "experiments.rds"))
}

#' @export
save_goldstandard = function(goldstandard, goldstandardid=goldstandard$info$id) {
  goldstandard_folder = file.path(.datasets_location, "goldstandards/", goldstandardid)
  dir.create(goldstandard_folder)
  saveRDS(goldstandard, file.path(goldstandard_folder, "goldstandard.rds"))
  
  goldstandards = readRDS(file.path(.datasets_location, "goldstandards.rds"))
  goldstandards %>% filter(id != goldstandardid) %>% bind_rows(goldstandard$info) %>% saveRDS(file.path(.datasets_location, "goldstandards.rds"))
}

#' @export
load_goldstandard = function(goldstandardid, contents=contents_goldstandard()) {
  goldstandard_folder = file.path(.datasets_location, "goldstandards/", goldstandardid)
  readRDS(file.path(goldstandard_folder, "goldstandard.rds"))
}
#' @export
delete_goldstandards = function(goldstandardids) {
  purrr::walk(goldstandardids, ~unlink(file.path(.datasets_location, "goldstandards/", .), T))
  goldstandards = readRDS(file.path(.datasets_location, "goldstandards.rds"))
  goldstandards %>% filter(!(id %in% goldstandardids))  %>% saveRDS(file.path(.datasets_location, "goldstandards.rds"))
}

#' @export
contents_dataset = function(counts=TRUE, platform=TRUE, experiment=contents_experiment(), goldstandard=contents_goldstandard(), model=contents_model(), info=TRUE) {
  as.list(environment(), all=TRUE)
}
#' @export
contents_goldstandard = function(goldstandard = TRUE) {
  as.list(environment(), all=TRUE)
}
#' @export
contents_model = function(model = TRUE) {
  as.list(environment(), all=TRUE)
}
#' @export
save_dataset = function(dataset, datasetid=dataset$info$id) {
  dataset_folder = file.path(.datasets_location, "datasets/", datasetid)
  dir.create(dataset_folder)
  purrr::walk(names(dataset), ~saveRDS(dataset[[.]], file.path(dataset_folder, paste0(., ".rds"))))
  
  datasets = readRDS(file.path(.datasets_location, "datasets.rds"))
  datasets %>% filter(id != datasetid) %>% bind_rows(dataset$info) %>% saveRDS(file.path(.datasets_location, "datasets.rds"))
}
#' @export
load_dataset = function(datasetid, contents = contents_dataset()) {
  combinedinfo = tibble(id=datasetid) %>% left_join(readRDS(file.path(".datasets_location", "datasets.rds")), by="id") %>% left_join(readRDS(file.path(".datasets_location", "goldstandards.rds")) %>% select(-version) %>% rename(goldstandardid=id), by="experimentid") %>% left_join(readRDS(file.path(".datasets_location", "experiments.rds")) %>% select(-version) %>% rename(experimentid=id), by="experimentid") %>% as.list()
  
  dataset = list()
  dataset_folder = file.path(.datasets_location, "/datasets/", datasetid)
  contents2 = keep(contents, ~!is_logical(.)) # all additional contents, not directly loaded from a dataset
  contents = keep(contents, ~is_logical(.)) # all contents loaded from a dataset
  loadin = names(contents)[purrr::map_lgl(contents, ~.)]
  dataset = purrr::map(loadin, ~readRDS(file.path(dataset_folder, paste0(., ".rds"))))
  names(dataset) = loadin
  dataset$experiment = load_experiment(dataset$info$experimentid, contents2$experiment)
  dataset$gs = load_goldstandard(combinedinfo$goldstandardid, contents2$goldstandard)
  dataset$model = load_model(combinedinfo$modelid, contents2$model)
  dataset
}


#' @export
refresh_experiments = function() readRDS(file.path(.datasets_location, "experiments.rds")) %>% filter(id %in% list.dirs(file.path(.datasets_location, "experiments/."), full.names=F)) %>% saveRDS(file.path(.datasets_location, "experiments.rds"))

#' @export
refresh_models = function() readRDS(file.path(.datasets_location, "models.rds")) %>% filter(id %in% list.dirs(file.path(.datasets_location, "models/."), full.names=F)) %>% saveRDS(file.path(.datasets_location, "models.rds"))

#' @export
refresh_datasets = function() readRDS(file.path(.datasets_location, "datasets.rds")) %>% filter(id %in% list.dirs(file.path(.datasets_location, "datasets/."), full.names=F)) %>% saveRDS(file.path(.datasets_location, "datasets.rds"))

#' @export
refresh_goldstandards = function() readRDS(file.path(.datasets_location, "goldstandards.rds")) %>% filter(id %in% list.dirs(file.path(.datasets_location, "goldstandards/."), full.names=F)) %>% saveRDS(file.path(.datasets_location, "goldstandards.rds"))

remove_duplicate_goldstandards = function() {
  readRDS(file.path(.datasets_location, "goldstandards.rds")) %>% group_by(experimentid) %>% filter(row_number() < n()) %>% .$id %>% delete_goldstandards()
}
#walk(c("experiments.rds", "models.rds", "datasets.rds", "goldstandards.rds"), function(x) {tibble() %>% write_rds(file.path(.datasets_location, x))})
