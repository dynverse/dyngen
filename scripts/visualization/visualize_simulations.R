datasets_info <- list_datasets()
datasets_info <- datasets_info %>% filter(stringr::str_detect(id, "^2017_04_25.*"))

experimentid = unique(datasets_info$experimentid)[[4]]


experiment = load_experiment(experimentid, contents_experiment(simulations=TRUE))
simulation_expressions = map(experiment$simulations, ~.$expression[-c(1:30), ])
cellinfo = tibble(simulationid=unlist(map(seq_along(simulation_expressions), ~rep(., nrow(simulation_expressions[[.]])))))
simulation_expressions = invoke(rbind, simulation_expressions)
simulation_expressions_modules = get_module_counts(simulation_expressions, experiment$model$modulemembership)


space = SCORPIUS::reduce.dimensionality.landmarked(simulation_expressions, SCORPIUS::correlation.distance, landmark.method = "naive", num.landmarks = 20)$S
space = space %>% as_tibble() %>% set_colnames(c("Comp1", "Comp2")) %>% bind_cols(cellinfo) %>% bind_cols(simulation_expressions_modules %>% as_tibble())

ggplot(space) + geom_point(aes(Comp1, Comp2))

plots = list()
plots = c(plots, list(ggplot(space) + geom_path(aes_string("Comp1", "Comp2", group="simulationid", color="factor(simulationid)")) + theme_nothing() + theme(legend.position="none") + scale_color_brewer(palette="Paired")))
plot_module = function(moduleid) {
  ggplot(space) + geom_path(aes_string("Comp1", "Comp2", group="simulationid", color=moduleid)) + viridis::scale_color_viridis(option="A", begin=0.1) + theme_nothing() + theme(legend.position="none")
}
plots_modules = colnames(simulation_expressions_modules) %>% map(plot_module)
cowplot::plot_grid(plotlist=plots_modules)
dev.off()
plots = c(plots, plots_modules[c(1, 4, 5, 10)])
#plots = c(plots, plots_modules[c(2, 3, 4, 5)])

save_plot("bifurcation_convergence.png", cowplot::plot_grid(plotlist=plots, ncol = 5), base_width = 12, base_height=3)
dev.off()


experiment$simulations