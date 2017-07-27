modulenetnames = c("linear", "bifurcating", "consecutive_bifurcating", "cycle", "trifurcating", "bifurcating_convergence", "bifurcating_loop", "bifurcating_cycle")

model = generate_model(modulenetname = "trifurcating")
experiment = run_experiment(model, 5, nsimulations=40)
experiment = run_experiment(model, 5, nsimulations=100)
experiment = run_experiment(model, 20, nsimulations=40)


gs = extract_goldstandard(experiment)
experiment = experiments[[5]]


gs$cellinfo = simulationstepinfo


filter_expression = function(expression, ncells=min(nrow(expression), 1000), Goi=colnames(expression)) {
  if(!is.null(ncells)) {
    expression = expression[sample(seq_len(nrow(expression)), ncells), ]
  }
  expression[, colnames(expression) %in% as.character(Goi)]
}



expression = experiment$expression %>% {log2(.+1)}
rownames(expression) = experiment$cellinfo$simulationstepid
source("~/thesis/projects/dyngen/scripts/evaluation/methods/dimred.R")
space = ica(expression, ndim = 2)
space = tsne(expression, ndim = 2)
space = mds(expression, ndim = 2)



plotdata = bind_cols(space %>% as.data.frame, simulationstepinfo %>% slice(match(rownames(space), simulationstepid)))
ggplot(plotdata) + geom_point(aes(Comp2, Comp1, color=factor(piecestateid)))# + viridis::scale_color_viridis(option="A")






combine_simulations <- function(simulations) {
  simulations <- smoothe_simulations(experiment$simulations, model)
  simulation_expressions = map(simulations, ~.$expression_smooth[seq(1, nrow(.$expression_smooth), 5), ])# %>% do.call(rbind, .)
  combined <- simulation_expressions %>% do.call(rbind, .)
  rownames(combined) = seq_len(nrow(combined))
  combined <- combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
  combinedinfo <- tibble(simulationid = factor(unlist(map(seq_along(simulation_expressions), ~rep(., times=nrow(simulation_expressions[[.]]))))), sample=seq_len(nrow(combined)), simulationstep = unlist(map(simulation_expressions, ~as.numeric(rownames(.)))))
  
  lst(combined, combinedinfo)
}
list2env(combine_simulations(experiment$simulations), .GlobalEnv)

space = mds_withlandmarks(combined, SCORPIUS::correlation.distance, landmark.method = "naive", k=2)$S %>% set_colnames(c("Comp1", "Comp2"))
space = ica(combined, ndim = 2)
plotdata = bind_cols(space %>% as.data.frame, combinedinfo)
ggplot(plotdata) + geom_path(aes(Comp2, Comp1, group=simulationid, color=simulationstep)) + viridis::scale_color_viridis(option="A")
ggplot(plotdata) + geom_path(aes(Comp2, Comp1, group=simulationid, color=simulationid))
