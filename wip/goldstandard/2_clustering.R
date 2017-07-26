simulations <- smoothe_simulations(experiment$simulations, model)


# first combine the simulations in one matrix (= combined) together with extra info about each "cell" (=combinedinfo)
simulation_expressions = map(simulations, ~.$expression_smooth[seq(1, nrow(.$expression_smooth), 5), ])# %>% do.call(rbind, .)
combined <- simulation_expressions %>% do.call(rbind, .)
rownames(combined) = seq_len(nrow(combined))
combined <- combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
combinedinfo <- tibble(simulationid = factor(unlist(map(seq_along(simulation_expressions), ~rep(., times=nrow(simulation_expressions[[.]]))))), sample=seq_len(nrow(combined)), simulationstep = unlist(map(simulation_expressions, ~as.numeric(rownames(.)))))

# check for convergence, important to check for cycles
check_convergence <- function(expression) {
  differences <- expression %>% diff(1) %>% abs %>% apply(1, mean)
  baseline <- mean(differences[seq(nrow(expression)-40, nrow(expression)-1)])
  
  (seq_len(nrow(expression)) >= (last(which(differences > baseline*5)) + 1)) %>% {ifelse(is.na(.), FALSE, .)}
}

combinedinfo$converged = unlist(simulation_expressions %>% map(check_convergence))

# calculate distances between each "cell"
distances <- as.matrix(proxy::dist(combined, method="euclidean"))
pheatmap(distances[1:400, 1:400], cluster_rows=F, cluster_cols=F)
distancesdf <- distances %>% reshape2::melt(value.name = "distance", varnames=c("from", "to")) %>% mutate(from_simulation=combinedinfo$simulationid[from], to_simulation = combinedinfo$simulationid[to])
distancesdf$to_simulation_edited = ifelse((distancesdf$from_simulation == distancesdf$to_simulation) & (abs(combinedinfo$simulationstep[distancesdf$from] - combinedinfo$simulationstep[distancesdf$to]) > 500), 0, distancesdf$to_simulation) # when there are enough steps between two "cells" of the same simulation, we should call this a "new" simulation as to detect cycles, here we assign to these cells a 0

cooccurences = distancesdf %>% reshape2::acast(from_simulation~to_simulation_edited, mean, value.var="distance")
cooccurences %>% pheatmap()

cooccurences = distancesdf %>% reshape2::acast(from~to_simulation_edited, mean, value.var="distance")
cooccurences %>% pheatmap()
cooccurences = distancesdf %>% reshape2::acast(from~to_simulation_edited, function(x) as.numeric(sum(x<1)), value.var="distance")# %>% {log2(.+1)}
#cooccurences = distancesdf %>% reshape2::acast(from~to_simulation_edited, function(x) as.numeric(sum(x<1)), value.var="distance") > 0


clust = hclust(proxy::dist(cooccurences, method="euclidean"), "ward.D2")
plot(clust)
labels = cutree(clust, k = 3)

labels = cluster::pam(proxy::dist(cooccurences, method="euclidean"), 5)$clustering

labels = MeanShift::bmsClustering(t(cooccurences), 40)$labels

filter_labels <- function(labels, minsize=50) {
  tab <- table(labels)
  wronglabels <- as.numeric(names(tab))[tab < minsize]
  labels[labels %in% wronglabels] <- NA
  as.numeric(factor(labels)) # reset ordering to 1, 2, 3, ... + NA
}

combinedinfo$bundle = filter_labels(labels) %>% zoo::na.locf(na.rm=FALSE)
ggplot(combinedinfo) + geom_line(aes(simulationstep, bundle, group=simulationid, color=simulationid))

bundlechanges = combinedinfo %>% {(.[seq_len(nrow(.)-1), ]$simulationid == .[seq_len(nrow(.)-1)+1, ]$simulationid) & (.[seq_len(nrow(.)-1), ]$bundle != .[seq_len(nrow(.)-1)+1, ]$bundle)}

statenet = tibble(from=combinedinfo$bundle[bundlechanges], to=combinedinfo$bundle[which(bundlechanges)+1]) %>% group_by(from, to) %>% summarise(nsimulations=n()) %>% ungroup()
statenet %>% graph_from_data_frame() %>% plot

# for cycles, check if this bundle reconnects with itself (ie. with the same simulation id but with multiple steps in between)
cycle_cutoff = 1000
statenet <- map(unique(combinedinfo$bundle), function(bundleid) {
  cycle_score = (cooccurences[(combinedinfo$bundle == bundleid) & (!combinedinfo$converged), "0"] %>% sum())
  print(cycle_score)
  if(cycle_score > cycle_cutoff) {
    tibble(from=bundleid, to=bundleid)
  }else {
    tibble()
  }
}) %>% bind_rows(statenet)

statenet %>% graph_from_data_frame() %>% plot


# labels2clusters <- function(labels, minsize=100) {
#   tab <- table(labels)
#   map(as.numeric(names(tab)[tab >= minsize]), ~as.numeric(names(labels)[labels == .]))
# }
# bundles = labels2clusters(labels)
# 
# bundles



library(dtw)














expression <- combined[combinedinfo$simulationid == 7, ]
distances <- as.matrix(proxy::dist(expression, combined))
#pheatmap(distances, cluster_rows=F, cluster_cols=F)

distance_cutoff = 2

repsimids = rep(combinedinfo$simulationid, 10)
prevpercs <- repsimids[(c(distances[1:10, ]) < distance_cutoff)] %>% table %>% {.>1}
splits = c()
percs = list()
for (i in seq_len(nrow(distances)-1-20)+1) {
  curpercs = repsimids[(distances[i:i+10, ] < distance_cutoff)] %>% table %>% {.>1}
  difference = (sum(abs(prevpercs - curpercs)))
  
  print(difference)
  
  percs = c(percs, list(curpercs))
  if(difference > 3) {
    prevpercs = repsimids[(distances[(i+10):(i+20), ] < distance_cutoff)] %>% table %>% {.>1}
    splits = c(splits, i)
  }
}
percs %>% map(~as.numeric(.)) %>% do.call(rbind, .) %>% pheatmap(cluster_rows=F)



differences = map(seq_len(nrow(distances)),  function(i) {
  combinedinfo[(distances[i, ] < 2), ]$simulationid %>% table %>% {.>1} %>% as.numeric()
}) %>% do.call(rbind, .)
differences %>% pheatmap(cluster_cols=T, cluster_rows=F, gaps_row=splits)


differences %>% pheatmap(cluster_cols=T, cluster_rows=F)
differences %>% as.matrix %>% diff(1) %>% abs %>% apply(1, sum) %>% plot



##





















## KNN is not the best, because what if simulations split up? Suddenly, the knn will look further away... Better to look at the distance in expression!
knn <- FNN::get.knn(combined, k = 100)
knn <- knn$nn.index %>% reshape2::melt(value.name = "to", varnames=c("from", "k"))
knn %>% mutate(from_simulation=combinedinfo$simulationid[from], to_simulation = combinedinfo$simulationid[to])

cooccurences = knn %>% mutate(from_simulation=combinedinfo$simulationid[from], to_simulation = combinedinfo$simulationid[to]) %>% reshape2::acast(from_simulation~to_simulation, length, value.var="to_simulation")
cooccurences %>% pheatmap()
