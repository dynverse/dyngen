


# now that we know to which 

get_mean_warp_times <- function(filtered) {
  pairs = expand.grid(a=seq_along(filtered$expression), b=seq_along(filtered$expression)) %>% as.data.frame() %>% as.data.frame() %$% map2(a, b, function(a, b) {
    tibble(a=a, b=b, dtw=list(dtw::dtw(filtered$expression[[a]], filtered$expression[[b]])))
  }) %>% bind_rows()
  
  match_to_other <- function(dt) {
    (dt$index2[match(unique(dt$index1), dt$index1)]/max(dt$index2))
  }
  
  times = map(seq_along(filtered$expression), function(simulationid) {
    pairs %>% filter(a==simulationid) %>% pull(dtw) %>% map(match_to_other) %>% invoke(rbind, .) %>% apply(2, mean) %>% set_names(rownames(filtered$expression[[simulationid]]))
  }) %>% unlist()
  
  times * mean(filtered$expression %>% map_int(~nrow(.)))
}

simulationstepinfo = pieces %>% split(pieces$piecestateid) %>% mclapply(get_mean_warp_times, mc.cores = 8) %>% map2(unique(pieces$piecestateid), ~tibble(time=.x, piecestateid=.y, simulationstepid=names(.x))) %>% bind_rows() %>% mutate(piecestateid=factor(piecestateid))

simulationstepinfo = pieces %>% mutate(simulationstepid = map(expression, ~rownames(.)), stepid = map2(start_stepid, end_stepid, ~seq(.x, .y)), piecestepid = map2(start_stepid, end_stepid, ~seq_len(.y-.x+1))) %>% unnest(simulationstepid, stepid, piecestepid) %>% group_by(pieceid) %>% mutate(time=piecestepid)


combined = map(experiment$simulations, ~.$expression) %>% do.call(rbind, .)
combined = combined[simulationstepinfo$simulationstepid, ]
space = ica(combined, ndim = 2)
space = tsne(combined, ndim = 2)
space = SCORPIUS::reduce.dimensionality.landmarked(combined, SCORPIUS::correlation.distance, landmark.method = "naive", k = 2, num.landmarks = 200)$S

plotdata = space %>% as.data.frame() %>% set_colnames(c("Comp1", "Comp2")) %>% bind_cols(simulationstepinfo) %>% filter(simulationid == 1)
ggplot(plotdata %>% mutate(piecestateid=factor(piecestateid))) + geom_point(aes(Comp2, Comp1, color=piecestateid, group=simulationid))
ggplot(plotdata) + geom_point(aes(Comp2, Comp1, color=time, group=simulationid), data=plotdata) + viridis::scale_color_viridis(option="A")

pheatmap(combined, cluster_cols=F, cluster_rows=F, annotation_row = plotdata %>% select(time))









smoothe_gaussian <- function(combined, times, timesoi=seq(0, 1, 0.02), sd=0.1) {
  map(timesoi, function(time) {
    densities = dnorm(times, time, sd)
    x = ((combined * densities) %>% colSums()) / sum(densities)
  }) %>% invoke(rbind, .) %>% {set_rownames(., seq_len(nrow(.)))} %>% list(expression=., cellinfo=data.frame(time=timesoi))
}

smoothed = smoothe_gaussian(combined, times)
pheatmap(smoothed$expression, cluster_cols=F, cluster_rows=F, annotation_row = smoothed$cellinfo)
