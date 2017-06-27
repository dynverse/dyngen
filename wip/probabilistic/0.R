endmodules <- c("M10", "M8", "M9")
endmodules <- c("M3", "M4")## consecutive bifurcating
#endmodules <- c("M4", "M5") ## bifurcating convergence

modulecounts <- lapply(seq_along(experiment$simulations), function(simulationid) {
  simulation <- experiment$simulations[[simulationid]]
  modulecounts <- get_module_counts(simulation$expression, model$modulemembership)
  
  simulationid %>% print
  ((abs(modulecounts[, endmodules]) %>% apply(., 1, max))/(abs(modulecounts[, endmodules]) %>% apply(., 1, min))) %>% plot
  
  ends <- ((abs(modulecounts[, endmodules]) %>% apply(., 1, max))/(abs(modulecounts[, endmodules]) %>% apply(., 1, min))) %>% {.>100} %>% which()
  if(length(ends) == 0) {
    return(data.frame())
  }
  
  firstend <- ends %>% .[[1]]
  endmodule <- endmodules[which.max(modulecounts[firstend, endmodules])]
  
  #endmodule <- endmodules[which.max(modulecounts[nrow(modulecounts), endmodules])]
  
  as.data.frame(modulecounts) %>% mutate(simulation=simulationid, endmodule=endmodule)
}) %>% bind_rows() %>% mutate(simulation=factor(simulation))

modulesoi <- list("M3"=c("M8", "M9"), "M4"=c("M10"))
modulesoi <- list("M7"=c("M8"), "M6"=c("M9"))
modulesoi <- list("M3"="M3", "M4"="M4") ## consecutive bifurcating
#modulesoi <- list("M4"="M4", "M5"="M5") ## bifurcating convergence

plotdata = experiment$expression_modules[,names(modulesoi)] %>% as.data.frame() %>% rownames_to_column("cell") %>% left_join(gs$cellinfo, by="cell")
ggplot(as.data.frame(modulecounts), aes_string(names(modulesoi)[[1]], names(modulesoi)[[2]])) + geom_path(aes(group=simulation, color=endmodule)) + coord_equal() + geom_point(aes(shape=piecestateid), data=plotdata)

filtermodulecounts1 = modulecounts %>% filter(endmodule %in% modulesoi[[1]])
filtermodulecounts2 = modulecounts %>% filter(endmodule %in% modulesoi[[2]])

endmodulecounts = modulecounts %>% group_by(simulation) %>% summarise(endmodule=first(endmodule)) %>% .$endmodule %>% table()

bw = 0.5
n = 100
lims = c(0,5,0,5)
densities1 = MASS::kde2d(filtermodulecounts1[, names(modulesoi)[[1]]], filtermodulecounts1[, names(modulesoi)[[2]]], bw, n, lims=lims)
xs = densities1$x
ys = densities1$y
densities1 = densities1$z# / endmodulecounts[[modulesoi[[1]]]]
#densities1 = densities1/sum(densities1)
densities2 = MASS::kde2d(filtermodulecounts2[, names(modulesoi)[[1]]], filtermodulecounts2[, names(modulesoi)[[2]]], bw, n, lims=lims)$z# / endmodulecounts[[modulesoi[[2]]]]
#densities2 = densities1/sum(densities2)
densities = MASS::kde2d(modulecounts[, names(modulesoi)[[1]]], modulecounts[, names(modulesoi)[[2]]], bw, n, lims=lims)$z
#pheatmap(log(densities+1), cluster_rows=F, cluster_cols=F)

densities_differences = (-densities2+densities1)/(densities1+densities2+0.0001)
pheatmap(densities_differences, cluster_rows=F, cluster_cols=F, color=colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100))
densities_differences_df = densities_differences %>% reshape2::melt(varnames=c("module1", "module2"), value.name="density_difference") %>% mutate(module1=xs[module1], module2=ys[module2])
densities_differences_df$density = as.numeric(densities)
plotdata = experiment$expression_modules[,names(modulesoi)] %>% as.data.frame() %>% rename_(module1=names(modulesoi)[[1]], module2=names(modulesoi)[[2]]) %>% rownames_to_column("cell") %>% left_join(gs$cellinfo, by="cell")
ggplot(densities_differences_df) + geom_raster(aes(module1, module2, fill=density_difference)) + scale_fill_gradientn(colors=c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'), breaks=seq(-1, 1, length.out=11)) + coord_equal() +
  geom_point(aes(module1, module2, color=progression), data=plotdata %>% filter(piecestateid == 1))


probable_cutoff <- 0.001
types <- tibble(type=c("notprobable", "notprobable", "uncomitted", "committed"), probable=c(F, F, T, T), comitted=c(F, T, F, T))
densities_differences_df <- densities_differences_df %>% mutate(
  probable = density > probable_cutoff, 
  comitted = abs(density_difference) > 0.9
) %>% left_join(types, by=c("probable", "comitted"))
ggplot(densities_differences_df) + geom_raster(aes(module1, module2, fill=type)) + coord_equal()




cellinfooi = gs$cellinfo %>% filter(piecestateid == 1)
cellinfooi$densitydifference = map_dbl(cellinfooi %>% .$cell, function(cell){
  row = experiment$expression_modules[cell,]
  densities_differences[
    which.min(abs(xs-row[names(modulesoi)[[1]]])),
    which.min(abs(xs-row[names(modulesoi)[[2]]]))
    ]
})


plot(cellinfooi$progression, cellinfooi$densitydifference)


cellinfooi$percA = 1- ((cellinfooi$progression) / max(gs$reference$cellinfo %>% filter(piecestateid == 1) %>% .$progression))
cellinfooi = cellinfooi %>% mutate(percB = (densitydifference + 1)/2 * (1-percA), percC = (1-percA) - percB)
plot(cellinfooi$percB, cellinfooi$percC)
cellinfooi = experiment$cellinfo %>% right_join(cellinfooi, by="cell")

g1 <- ggplot(cellinfooi) + geom_point(aes(percB, percC, color=percA)) + viridis::scale_color_viridis(option="A") + coord_equal()
g1

ggplot(experiment$cellinfo %>% right_join(cellinfooi, by="cell") %>% arrange(simulationtime)) + geom_path(aes(percA, percB-percC, group = simulationid)) + viridis::scale_color_viridis(option="A") + coord_equal()

ggplot(cellinfooi %>% arrange(simulationtime)) + geom_path(aes(percA, percB-percC, colour = factor(simulationid))) + coord_equal()

ggplot(cellinfooi %>% arrange(simulationtime) %>% mutate(gr = floor(simulationid / 5), gr2 = simulationid %% 5)) + geom_path(aes(percA, percB-percC, colour = factor(gr2))) + coord_equal() + facet_wrap(~gr)

