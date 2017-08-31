simulations <- smoothe_simulations(experiment$simulations, model)
tobecomes <- extract_tobecomes(simulations, model$piecestates)
pieces <- divide_pieces(simulations, tobecomes, model$piecestates)

simulationstepinfo <- pieces %>% unnest(cellid, stepid, piecestepid) %>% select(-start_stepid, -end_stepid)

simulationstepinfo %<>% mutate(time = piecestepid)

##
combined = map(experiment$simulations, ~.$expression) %>% do.call(rbind, .)
expression = combined[simulationstepinfo$cellid, ]

expression = experiment$expression %>% set_rownames(experiment$cellinfo$simulationstepid)
##

space = ica(expression, ndim = 2)
space = SCORPIUS::reduce.dimensionality.landmarked(expression, SCORPIUS::correlation.distance, landmark.method = "naive", k = 2, num.landmarks = 200)$S


plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% mutate(cellid=rownames(.))
ggplot(plotdata) + geom_point(aes(Comp2, Comp1))





plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% mutate(cellid=rownames(.)) %>% left_join(gs$cellinfo, by="cellid")
ggplot(plotdata %>% mutate(piecestateid=factor(piecestateid))) + geom_point(aes(Comp2, Comp1, color=piecestateid, group=simulationid))
