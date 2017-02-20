imagefolder = "/home/wouters/thesis/presentations/1702_dambi/img/"

# Three problems
# 1: Start and end state of circular trajectory
# 2: Missing tail trajectory of dying cells
# 3: Pseudo time warping


celltypeoi = 2
Ec1 =  Emrna2[paste0("x_", genes %>% filter(celltype == 1) %>% .$gene), exprs(E)[paste0("rg_", 1),] > 0.2]
Ec2 =  Emrna2[paste0("x_", genes %>% filter(celltype == 2) %>% .$gene), exprs(E)[paste0("rg_", 2),] > 0.2]
Ecs = list(Ec1, Ec2)

Ecs = map(Ecs, ~.[apply(exprs(.), 1, sd) > 0,apply(exprs(.), 2, sd) > 0])
save(Ecs, file="Ecs.RData")
load(file="Ecs.RData")

Ec1 = Ecs[[1]]
Ec2 = Ecs[[2]]

pheatmap(SCORPIUS::quant.scale(exprs(Ec2[,phenoData(Ec2)$time %>% order]) %>% t, 0.05) %>% t, cluster_cols = F, scale="none", cluster_row=T, show_rownames = F, show_colnames = F)

spaces = map(Ecs, ~SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(t(exprs(.))),ndim = 3))
trajectories = map(spaces, ~SCORPIUS::infer.trajectory(.))

set.seed(4)
map(1:2, function(i) {
  plotdata = as.data.frame(spaces[[i]])
  plotdata$progression = phenoData(Ecs[[i]])$time
  ggplot(plotdata) + geom_jitter(aes(Comp1, Comp2, color=progression), width=.01, height=.01) + viridis::scale_color_viridis(end=0.8, option = "A", direction = -1) + theme_classic() + geom_path(aes(Comp1, Comp2), data=trajectories[[i]]$path %>% as.data.frame)
}) %>% cowplot::plot_grid(plotlist=.)# %>% cowplot::save_plot(file.path(imagefolder, "interacting_trajectories.png"), ., base_width=12, base_height=6)

trajectories[[1]] = list(time=spaces[[1]][,1] %>% {(. - min(.)) / (max(.) - min(.))})


rgl::plot3d(spaces[[2]], col=viridis::magma(20)[phenoData(Ecs[[2]])$time %>% cut(20, labels=F)])

pheatmap(SCORPIUS::quant.scale(t(exprs(Ec1)[,phenoData(Ec1)$time %>% order]), 0.05) %>% t, cluster_cols = F, scale="none", cluster_rows=F)



## make bulk expression

generate_bulk = function(Ec1, Ec2, sd=0.1, totaltimeoi=8, ntimepoints=10, verbose=F) {
  timepoints = seq(0, 1, length.out=ntimepoints)[-ntimepoints]
  Ec2_bulk = sapply(timepoints, function(timepoint) {
    weights = dnorm(phenoData(Ec2)$time / totaltimeoi, timepoint, sd)
    apply(exprs(Ec2), 1, function(row) weighted.mean(row, weights))
  })
  # Ec2_bulk %>% pheatmap(scale="row", cluster_cols=F)
  
  Ec1_bulk = sapply(timepoints, function(timepoint) {
    weights = dnorm(phenoData(Ec1)$time / totaltimeoi, timepoint, sd)
    if(sum(weights) < ncol(Ec1)/2) {  ## remove if not enough cells
      return(runif(nrow(Ec1), 0, 0.0001))
    }
    apply(exprs(Ec1), 1, function(row) weighted.mean(row, weights))
  })
  # Ec1_bulk %>% pheatmap(scale="row", cluster_cols=F)
  
  list(Ec1_bulk=Ec1_bulk, Ec2_bulk=Ec2_bulk)
}
sd = 0.05
timepoints = seq(0, 1, length.out=10)[-10]
lapply(timepoints, function(timepoint) tibble(timepoint=timepoint, time=seq(0, 1, 0.01), density=dnorm(seq(0, 1, 0.01), timepoint, sd))) %>% bind_rows %>% mutate(timepoint=factor(timepoint)) %>% ggplot() + geom_ribbon(aes(time, ymax=density, group=timepoint, color=timepoint, fill=timepoint), ymin=0, alpha=0.1)
lapply(timepoints, function(timepoint) tibble(timepoint=timepoint, time=seq(0, 2*pi, 0.01), density=circular::dwrappednormal(seq(0, 2*pi, 0.01), timepoint*2*pi, sd=sd*2*pi))) %>% bind_rows %>% mutate(timepoint=factor(timepoint), time=time/2/pi) %>% ggplot() + geom_ribbon(aes(time, ymax=density, group=timepoint, color=timepoint, fill=timepoint), ymin=0, alpha=0.1)


Ec_bulk = generate_bulk(Ec1, Ec2)
Ec1_bulk = Ec_bulk$Ec1_bulk
Ec2_bulk = Ec_bulk$Ec2_bulk

phenoData(Emrna)$time %>% hist


map_common_time = function(Ec1, Ec1_bulk, Ec2, Ec2_bulk) {
  ## smooth
  commonpseudotimes = seq(0, 1, length.out=200)
  Ec1_bulk_smoothed = apply(Ec1_bulk, 1, function(row) smooth.spline(seq(0, 1, length.out=length(row)), row) %>% predict(commonpseudotimes) %>% .$y) %>% t
  # Ec1_bulk_smoothed %>% pheatmap(cluster_cols=F, scale="row")
  
  Ec2_bulk_smoothed = apply(Ec2_bulk, 1, function(row) smooth.spline(seq(0, 1, length.out=length(row)), row) %>% predict(commonpseudotimes) %>% .$y) %>% t
  # Ec2_bulk_smoothed %>% pheatmap(cluster_cols=F, scale="row")
  
  ## actual mapping
  
  mapping1 = cor(exprs(Ec1), Ec1_bulk_smoothed) %>% apply(1, which.max)
  # plot(
  #   Ec1[,names(mapping1)] %>% phenoData() %>% .$time,
  #   commonpseudotimes[mapping1]
  # )
  # plot(
  #   trajectories[[1]]$time[names(mapping1)],
  #   commonpseudotimes[mapping1]
  # )
  
  mapping2 = cor(exprs(Ec2), Ec2_bulk_smoothed) %>% apply(1, which.max)
  # plot(
  #   Ec2[,names(mapping2)] %>% phenoData() %>% .$time,
  #   commonpseudotimes[mapping2]
  # )
  # plot(
  #   trajectories[[2]]$time[names(mapping2)],
  #   commonpseudotimes[mapping2]
  # )
  
  allcells = union(sampleNames(Ec1), sampleNames(Ec2))
  mappingdata1 = tibble(cell=names(mapping1))
  mappingdata1$original1 = trajectories[[1]]$time[mappingdata1$cell]
  mappingdata1$common1 = commonpseudotimes[mapping1]
  mappingdata1$simulation1 = (phenoData(Ec1)$time %>% setNames(sampleNames(Ec1)))[mappingdata1$cell]
  
  mappingdata2 = tibble(cell=names(mapping2))
  mappingdata2$original2 = trajectories[[2]]$time[mappingdata2$cell]
  mappingdata2$common2 = commonpseudotimes[mapping2]
  mappingdata2$simulation2 = (phenoData(Ec2)$time %>% setNames(sampleNames(Ec2)))[mappingdata2$cell]
  
  mappingdata = full_join(mappingdata1, mappingdata2, "cell")
  
  mappingdata
}

mappingdata = map_common_time(Ec1, Ec1_bulk, Ec2, Ec2_bulk)

#####

get_timediff = function(time1, time2, verbose=F) {
  differences = pmin((abs(time1 - time2) - 1) %>% abs, (time1 - time2) %>% abs)

  plotdata = tibble(time1, time2, cell=1:length(time1), difference=differences) %>% gather(variable, time, -cell, -difference)
  if(verbose) {
    plot(ggplot(plotdata) + 
      geom_line(aes(time, variable, group=cell, color=difference)) +
      geom_point(aes(time, variable)) +
      viridis::scale_color_viridis(limits=c(0, 0.5))
    )
  }
  
  differences %>% keep(!is.na(.))
}

get_timediff(mappingdata$simulation1, mappingdata$simulation2, T)
ggsave(file.path(imagefolder, "mapping_simulation.png"))
get_timediff(mappingdata$original1, mappingdata$original2, T)
ggsave(file.path(imagefolder, "mapping_original.png"))
get_timediff(mappingdata$common1, mappingdata$common2, T)


score_timediff = function(time1, time2, verbose=F) {
  differences=get_timediff(time1, time2, verbose)
  
  if (verbose) differences %>% hist
  
  differences %>% mean
}
score_timediff(mappingdata$original1, mappingdata$original2)
score_timediff(mappingdata$common2, mappingdata$common1)


score_timediff_variance = function(time1, time2, verbose=F) {
  differences=get_timediff(time1, time2, verbose)
  sd(differences)
}

score_timediff_variance(mappingdata$simulation1, mappingdata$simulation2)
score_timediff_variance(mappingdata$original1, mappingdata$original2)
score_timediff_variance(mappingdata$common2, mappingdata$common1)




####
allcells = union(sampleNames(Ec1), sampleNames(Ec2))
score = function(paramsetting) {
  if(paramsetting$timesource == "simulation") {
    time1 = phenoData(Ec1)$time %>% setNames(sampleNames(Ec1)) %>% .[allcells]
    time2 = phenoData(Ec2)$time %>% setNames(sampleNames(Ec2)) %>% .[allcells]
    c(list(
      timediff=score_timediff(time1, time2),
      timediff_variance = score_timediff_variance(time1, time2)
    ), paramsetting)
  } else if (paramsetting$timesource == "original") {
    time1 = trajectories[[1]]$time %>% setNames(sampleNames(Ec1)) %>% .[allcells]
    time2 = trajectories[[2]]$time %>% setNames(sampleNames(Ec2)) %>% .[allcells]
    c(list(
      timediff=score_timediff(time1, time2),
      timediff_variance = score_timediff_variance(time1, time2)
    ), paramsetting)
  } else {
    Ec_bulk = generate_bulk(Ec1, Ec2, paramsetting$sd, ntimepoints = paramsetting$ntimepoints)
    mappingdata = map_common_time(Ec1, Ec_bulk$Ec1_bulk, Ec2, Ec_bulk$Ec2_bulk)
    c(list(
      timediff=score_timediff(mappingdata$common2, mappingdata$common1),
      timediff_variance = score_timediff_variance(mappingdata$common2, mappingdata$common1)
    ), paramsetting)
  }
}


paramsettings = expand.grid(
  sd=c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
  ntimepoints=c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 100, 1000)
)
paramsettings$timesource = "common"
paramsettings = bind_rows(paramsettings,
  list(timesource="original"),
  list(timesource="simulation")
)


scores = mclapply(1:nrow(paramsettings), function(paramsettingid) {
  paramsetting = paramsettings[paramsettingid, ,drop=T]
  score(paramsetting)
}, mc.cores = 8) %>% bind_rows
scores %>% ggplot() + geom_raster(aes(factor(sd), factor(ntimepoints), fill=timediff))


plots = list()
bind_rows(scores %>% filter(sd == 0.05 & ntimepoints==8), scores %>% filter(timesource != "common")) %>% ggplot() + geom_boxplot(aes(timesource, timediff), data=scores %>% filter(timesource== "common"), alpha=0.1) + geom_point(aes(timesource, timediff))
plots = c(plots, list(last_plot() + ggtitle("Time difference")))

bind_rows(scores %>% filter(sd == 0.05 & ntimepoints==8), scores %>% filter(timesource != "common")) %>% ggplot() + geom_boxplot(aes(timesource, timediff_variance), data=scores %>% filter(timesource== "common"), alpha=0.1) + geom_point(aes(timesource, timediff_variance))
plots = c(plots, list(last_plot() + ggtitle("Time warping")))

cowplot::plot_grid(plotlist=plots, nrow=2) %>% cowplot::save_plot(file.path(imagefolder, "mapping_scoring.png"), ., base_width=4, base_height=8)




plots = list()
get_timediff(mappingdata$simulation1, mappingdata$simulation2, T)
plots = c(plots, list(last_plot() + ggtitle("Simulation time")))
get_timediff(mappingdata$original1, mappingdata$original2, T)
plots = c(plots, list(last_plot() + ggtitle("Original pseudo-time")))
get_timediff(mappingdata$common1, mappingdata$common2, T)
plots = c(plots, list(last_plot() + ggtitle("Common pseudo-time")))
cowplot::plot_grid(plotlist=plots, nrow=3) %>% cowplot::save_plot(file.path(imagefolder, "mappings.png"), ., base_width=8, base_height=8)
