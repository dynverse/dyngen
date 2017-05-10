experiments_info <- readRDS("results/experiments.rds")
experiments_info <- experiments_info %>% filter(stringr::str_detect(id, "^2017_04_25.*"))
.datasets_location <- "/home/wouters/thesis/projects/dyngen/results/"
experiments <- map(experiments_info$id, ~dyngen::load_experiment(., contents_experiment(T, T, T, T, T, T, T)))

experiment <- experiments[[1]] # select experiment

expressions <- list()
cellinfo <- tibble()
for (simulationid in seq_along(experiment$simulations)) {
  expression <- experiment$simulations[[simulationid]]$expression
  #expression <- expression[sort(sample(seq_len(nrow(expression)), 500)), ]
  rownames(expression) <- paste0(simulationid, "#", rownames(expression))
  cellinfo <- bind_rows(cellinfo, tibble(simulationid=simulationid, cell=rownames(expression)))
  expressions <- c(expressions, list(expression))
}
expression <- do.call(rbind, expressions)
cellinfo <- cellinfo %>% mutate(simulationid=factor(simulationid))

space <- dambiutils::mds_withlandmarks(expression, SCORPIUS::correlation.distance, num.landmarks=100, landmark.method="naive")$S
#space <- SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(expression))
space2 <- bind_cols(cellinfo, space %>% as.data.frame())
ggplot(space2) + geom_path(aes(V1, V2, group=simulationid, color=simulationid))
ggplot(space2) + geom_path(aes(V1, V2, group=simulationid, color=simulationid)) + facet_wrap(~simulationid)





modulecounts = dyngen::get_module_counts(expression, experiment$model$modulemembership) %>% reshape2::melt(varnames=c("cell", "module"), value.name="expression")
space3 <- left_join(space2, modulecounts, by="cell")
ggplot(space3) + geom_path(aes(V1, V2, group=simulationid, color=expression)) + facet_wrap(~module) + viridis::scale_color_viridis(option="A")
