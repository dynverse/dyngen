library(tidyverse)
library(dynutils)
#library(dyngen)

updates <- tribble(
  ~modulenet_name, ~totaltime,
  "linear", 5,
  "bifurcating", 5,
  "linear_long", 30,
  "cycle", 30,
  "consecutive_bifurcating", 10,
  "bifurcating_convergence", 15,
  "trifurcating", 10,
  "bifurcating_loop", 30
)



settings <- pmap(
  updates,
  function(modulenet_name, totaltime) {
    modifyList(dyngen:::simple_params, list(model = list(modulenet_name = modulenet_name), simulation=list(totaltime = totaltime)))
  }
)

ncores <- 8

params <- settings[[4]]

outputs <- PRISM::qsub_lapply(qsub_config=PRISM::override_qsub_config(num_cores = ncores, memory="4G"), qsub_environment=list2env(list()), settings, function(params) {
  library(tidyverse)
  
# outputs <- pbapply::pblapply(settings, function(params) {
  options(ncores = 8)
  
  # model
  model <- invoke(dyngen:::generate_model_from_modulenet, params$model)
  
  # simulation
  simulation <- invoke(dyngen:::simulate_multiple, params$simulation, model$system)
  
  simulation <- dyngen:::preprocess_simulation_for_gs(simulation, model, params$gs$smooth_window) # do preprocessing separate, otherwise zoo will stay in an infinite loop in case of later error
  gs <- invoke(dyngen:::extract_goldstandard, params$gs, simulation, model, preprocess=FALSE)
  # dyngen:::plot_goldstandard(simulation, model, gs)
  gs$checks <- dyngen:::check_goldstandard(gs)
  
  lst(simulation, model, gs)
})


list2env(outputs[[4]], .GlobalEnv)

pdf("hi.pdf", width=12, height=12)
walk(outputs, function(output) {
  list2env(output, environment())
  dyngen:::plot_net(model, label=FALSE, main_only = FALSE)
  dyngen:::plot_modulenet(model)
  dyngen:::plot_goldstandard(simulation, model, gs)
})
dev.off()



