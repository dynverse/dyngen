base_model = load_circuit("cell_cycle_arabidopsis")
models = map(1:1000, function(i) {
  model = base_model
  
  model = generate_formulae(model$net, model$geneinfo, model$celltypes) %>% c(model)
  model = generate_kinetics(model$vargroups, model$variables, model$nus.changes) %>% c(model) 
  model$burngenes = determine_burngenes(model)
  
  model
})


fitnesses = mclapply(models, function(model) {
  experiment = run_cycle(model)
  fitness_cycle(experiment)
}, mc.cores = 8) %>% unlist()

plot(fitnesses)

model = models[[which.max(fitnesses)]]
experiment = run_cycle(model)
