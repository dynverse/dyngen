library(tidyverse)
library(dynutils)
#library(dyngen)

updates <- tribble(
  ~modulenet_name, ~totaltime,
  "linear", 5,
  "bifurcating", 5,
  "linear_long", 20,
  "cycle", 30,
  "consecutive_bifurcating", 10,
  "bifurcating_convergence", 15,
  "trifurcating", 10,
  "bifurcating_loop", 30
)



settings <- pmap(
  updates,
  function(modulenet_name, totaltime) {
    modifyList(dyngen:::base_params, list(model = list(modulenet_name = modulenet_name), simulation=list(totaltime = totaltime)))
  }
)

params <- settings[[2]]

outputs <- map(settings, function(params) {
  params$simulation$local <- 8
  
  # model
  model <- invoke(dyngen:::generate_model_from_modulenet, params$model)
  dyngen:::plot_net(model, label=FALSE, main_only = FALSE)
  dyngen:::plot_modulenet(model)
  
  # simulation
  simulation <- invoke(simulate_multiple, params$simulation, model$system)
  
  list(simulation=simulation, model=model)
})
names(outputs) <- modulenet_names

list2env(outputs$consecutive_bifurcating, .GlobalEnv)

simulation <- preprocess_simulation_for_gs(simulation, model, params$gs$smooth_window) # do preprocessing separate, otherwise zoo will stay in an infinite loop in case of later error
gs <- invoke(extract_goldstandard, params$gs, simulation, model, preprocess=FALSE)
plot_goldstandard(simulation, model, gs)
check_goldstandard(gs)
