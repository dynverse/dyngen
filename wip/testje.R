modulenetname="linear"; totaltime=5; burntime=2; nsimulations = 40; ncells = 500
modulenetname="cycle"; totaltime=40; burntime=2; nsimulations = 40; ncells = 500
modulenetname="consecutive_bifurcating"; totaltime=8; burntime=2; nsimulations = 40; ncells = 500
modulenetname="bifurcating_convergence"; totaltime=8; burntime=2; nsimulations = 40; ncells = 500
modulenetname="trifurcating"; totaltime=5; burntime=2; nsimulations = 40; ncells = 500


experiment = run_experiment(modulenetname, totaltime, burntime, nsimulations, ncells)
dataset = run_scrnaseq(experiment, platforms[1, ] %>% as.list)
dataset$gs = extract_goldstandard(experiment, verbose=T)



plot(dataset$gs$cellinfo$progression, dataset$cellinfo$simulationtime)
plot_piecestate_changes(dataset)
plot_piecestate_progression(dataset)



scriptfolder = "./scripts/simulation/"
walk(c("6_goldstandard2.R"), ~source(file.path(scriptfolder, .)))

space = mds(dataset$expression) %>% as.data.frame
space$distance = dataset$gs$celldistances[501, ]
space = space %>% bind_cols(dataset$gs$cellinfo)

space %>% ggplot() + geom_point(aes(Comp1, Comp2, color=distance)) + viridis::scale_color_viridis()
