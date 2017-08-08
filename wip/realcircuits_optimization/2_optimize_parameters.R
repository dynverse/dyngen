library(GA)
library(tidyverse)
library(igraph)

base_model = generate_model("cycle")
base_model = load_circuit("cell_cycle_arabidopsis")

plot_net(base_model)

base_model = generate_formulae(base_model$net, base_model$geneinfo, base_model$celltypes) %>% c(base_model)
base_model = generate_kinetics(base_model$vargroups, base_model$variables, base_model$nus.changes) %>% c(base_model) 
base_model$burngenes = determine_burngenes(base_model)

name2gene = base_model$geneinfo$gene %>% set_names(base_model$geneinfo$name)

degradations = tibble()


varsettings = list(
  list(type="k", min=0.001, max=2),
  list(type="r", min=1, max=10),
  list(type="d", min=1, max=10),
  list(type="p", min=0.2, max=2),
  list(type="q", min=0.2, max=2),
  list(type="a0", min=0, max=1),
  list(type="a", min=0, max=1)
) %>% bind_rows()
#vargroupsoi = c("r", "d", "p", "q")
vargroupsoi = c("k", "a0", "a")
dynamic_params = unlist(base_model$vargroups[vargroupsoi])
static_params = setdiff(names(base_model$params), dynamic_params)
dynamic_paramsettings = tibble(
    type = unlist(map2(base_model$vargroups[vargroupsoi], vargroupsoi, ~rep(.y, length(.x)))), 
    original = base_model$params[dynamic_params]
  ) %>% left_join(varsettings, "type")
static_paramvalues = base_model$params[static_params]

fit = fit_cycle

newparams = dynamic_paramsettings$original
newparams = dynamic_paramsettings$original * 20
newparams = result@solution[1, ]
model1 = adapt_model(newparams, base_model)
experiment = run_cycle(model1)
fitness_cycle(experiment)

result = ga(type="real-valued", fitness=fit, min=dynamic_paramsettings$min, max=dynamic_paramsettings$max, maxiter=10, popSize=20, parallel=T, suggestions = dynamic_paramsettings$original)

result = qsub_run(
  function(x) {
    library(dyngen)
    GA::ga(type="real-valued", fitness=fit, min=dynamic_paramsettings$min, max=dynamic_paramsettings$max, maxiter=10, popSize=1000, parallel=T)
  }, 
  override_qsub_config(num_cores=24)
)[[1]]


experiment$expression %>% hist
space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(experiment$expression))
plotdata = space %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% left_join(as.data.frame(experiment$expression_modules) %>% tibble::rownames_to_column("cell"), by="cell") %>% bind_cols(experiment$cellinfo)
ggplot(plotdata) + geom_point(aes(Comp1, Comp2, color=simulationtime)) + viridis::scale_color_viridis()

pheatmap(experiment$expression, cluster_cols=F, cluster_rows=F)


model




expression_smooth = experiment$expression_modules %>% zoo::rollmean(50, c("extend", "extend", "extend")) %>% set_rownames(rownames(experiment$expression))
