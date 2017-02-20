library(GENIE3)

shift_expression = function(E, delay, get_mapping=F) {
  mapping = map_int(phenoData(E)$time - delay, function(time) which.min(abs(phenoData(E)$time - time)))
  if(!get_mapping) {
    E[,mapping]
  } else {
    mapping
  }
}

# ni
testbag = paste0("x_", sample(allgenes, length(allgenes)/2))

delay = 2
Eregulators = shift_expression(Emrna, delay)[,order(phenoData(Emrna)$time)]
Etargets = Emrna[,order(phenoData(Emrna)$time)]

result = genie3(t(exprs(Etargets)), regulators = paste0("x_", tfs), mc.cores = 8)
result2 = cor(t(exprs(Etargets))) %>% .[paste0("x_", tfs),] %>% abs
result3 = cor(t(exprs(Eregulators)), t(exprs(Etargets))) %>% .[paste0("x_", tfs),] %>% abs
source("scripts/1/lagged_genie3.R")
result4 = genie3(t(exprs(Etargets)), regulators = paste0("x_", tfs), data.regulators = t(exprs(Eregulators)), mc.cores = 8)
result5 = genie3(t(exprs(Etargets)), regulators = paste0("x_", tfs), data.regulators = t(exprs(Eregulators2)), mc.cores = 8)

# performance
cal_performance = function(observed_net, net, exclude=c(), include=union(colnames(observed_net), rownames(observed_net)), plot=F) {
  net$from = paste0("x_",net$from)
  net$to = paste0("x_",net$to)
  
  observed_net = observed_net[setdiff(intersect(rownames(observed_net), include), exclude), setdiff(intersect(colnames(observed_net), include), exclude)]
  net = net[(net$from %in% include) & !(net$from %in% exclude) & (net$to %in% include) & !(net$to %in% exclude),]
  net$from = factor(net$from, levels=rownames(observed_net))
  net$to = factor(net$to, levels=colnames(observed_net))
  
  net.adj = net %>% reshape2::acast(from~to, function(x) any(x>0), value.var="strength", drop=F)
  
  #net.adj = net.adj[rownames(observed_net), colnames(observed_net)]
  
  prediction = ROCR::prediction(observed_net %>% as.numeric, net.adj %>% as.logical)
  performance1 = ROCR::performance(prediction, measure="tpr", x.measure="fpr")
  performance2 = ROCR::performance(prediction, measure="prec", x.measure="rec")
  performance = tibble(tpr = performance1@y.values[[1]], fpr=performance1@x.values[[1]], precision=performance2@y.values[[1]], recall=performance2@x.values[[1]])
  
  performance = performance %>% replace_na(list(recall=0, precision=1))
  aupr = caTools::trapz(performance$recall, performance$precision)
  
  #performance %>% ggplot() + geom_line(aes(fpr, tpr))
  #performance %>% arrange(precision) %>% ggplot() + geom_point(aes(recall, precision))
  
  if (plot) limma::barcodeplot(result3 %>% as.numeric, net.adj)
  
  list(aupr=aupr, auc=ROCR::performance(prediction, measure="auc")@y.values[[1]])
}

bind_rows(lapply(list(result, result2, result3, result4, result5), function(result) cal_performance(result, net, exclude=testbag)))
### find the delay

find_delay = function(from, to, plot=F) {
  crosscor = ccf(exprs(E)[from,], exprs(E)[to,], lag.max = 100, plot=plot)
  crosscor$lag[,1,1][[which.max(crosscor$acf[,1,1])]]
}
find_delay = function(from, to, maxdelay = 100, plot=F) {
  E = Emrna[,sampleNames(Emrna)[order(phenoData(Emrna)$time)]]
  crosscor = lapply(seq(0, maxdelay), function(delay) cor(exprs(E)[from,1:(ncol(E)-delay)], exprs(E)[to,(delay+1):ncol(E)]))
  mean(diff(sort(phenoData(E)$time))) * (which.max(crosscor) - 1)
}


net = net %>% rowwise() %>% mutate(delay=find_delay(from, to, maxdelay = 50))
delay = mean(net$delay)
net %>% ggplot() + geom_histogram(aes(delay)) + geom_vline(xintercept=delay)

linkid = 2
from = paste0("x_", net[linkid,]$from)
to = paste0("x_", net[linkid,]$to)

plotdata = data.frame(t(exprs(Emrna)[c(from, to),]), time=phenoData(Emrna)$time) %>% rename_(from=from, to=to)
plotdata %>% 
  ggplot() + geom_smooth(aes(time, from), color="green", span=0.5) + geom_smooth(aes(time, to), color="blue", span=0.5) + geom_smooth(aes(time-delay, to), color="#232784", span=0.5)+ geom_point(aes(time, from), color="green", alpha=0.2) + geom_point(aes(time, to), color="blue", alpha=0.2) + geom_point(aes(time-delay, to), color="#232784", alpha=0.2)


delay = 2
cor(exprs(Emrna)[from,1:(ncol(Emrna)-delay)], exprs(Emrna)[to,(delay+1):ncol(Emrna)])
plot(exprs(Emrna)[from,1:(ncol(Emrna)-delay)], exprs(Emrna)[to,(delay+1):ncol(Emrna)])
plot(exprs(Emrna)[from,], exprs(Emrna)[to,])


Emrna_shifted = shift_expression(Emrna, 3)
plotdata = tibble(
  from_original=exprs(Emrna)[from,],
  to_original=exprs(Emrna)[to,],
  from_shifted=exprs(Emrna_shifted)[from,],
  to_shifted=exprs(Emrna)[to,]
)
plotdata %>% ggplot() + 
  geom_point(aes(from_original, to_original), color="blue") + 
  geom_point(aes(from_shifted, to_shifted), color="red")
plot(plotdata$from_shifted, plotdata$to_shifted)

######
allgenes = c(tfs, targets)
netmask = net %>% reshape2::acast(from~to, function(x) any(x>0), value.var="strength")

ncells = ncol(Emrna)
mapping = c(1:ncells)
curscore = 0

steps = list()

# method 1: in each iteration try move each cell upward and choose the most optimal one
for (i in 1:2) {
  newmapping_scores = lapply(1:(ncells-1), function(cellid) {
    if(mapping[[cellid]] == mapping[[cellid+1]]) return(list(auc=0, aupr=0, newmapping=list(mapping)))
    newmapping = mapping
    newmapping[[cellid]] = newmapping[[cellid]] + 1
    Eregulators = Etargets[,newmapping]
    
    observed_net = cor(t(exprs(Eregulators[paste0("x_", tfs),])), t(exprs(Etargets))) %>% abs
    
    scores = cal_performance(observed_net, net)
    scores$newmapping = list(newmapping)
    scores
  }) %>% bind_rows
  
  newmapping_id = newmapping_scores$auc %>% which.max
  newscore = newmapping_scores$auc[[newmapping_id]]
  print(newscore - curscore)
  if (newscore > curscore) {
    mapping = newmapping_scores$newmapping[[newmapping_id]]
  } else {
    break()
  }
  steps[[length(steps) + 1]] = list(original=1:length(mapping), mapping=mapping, stepid=i, auc=newscore)
}

# method 2: from the last cell, try to move each cell as far as possible as long as it increases performance
# in principle you could allow a cell to go further than the next cell, eg. to improve TI, although this is not advisable IMO
stepid = 0
maxstepsize = 50
for (i in 1:1) {
  for (cellid in 2:ncells) {
    newmapping_scores = lapply(0:maxstepsize, function(stepsize) {
      if (mapping[[cellid]]-stepsize < 1) return(list(auc=0, aupr=0, newmapping=list(mapping)))
      #if ((mapping[[cellid]]-stepsize < 1) | (mapping[[cellid-1]] > mapping[[cellid]] - stepsize)) return(list(auc=0, aupr=0, newmapping=list(mapping)))
      newmapping = mapping
      newmapping[[cellid]] = newmapping[[cellid]] - stepsize
      Eregulators = Etargets[paste0("x_", tfs), newmapping]
      
      observed_net = cor(t(exprs(Eregulators)), t(exprs(Etargets))) %>% abs
      
      scores = cal_performance(observed_net, net)
      scores$newmapping = list(newmapping)
      scores
    }) %>% bind_rows
    newmapping_id = newmapping_scores$auc %>% which.max
    newauc = newmapping_scores$auc[[newmapping_id]]
    newaupr = newmapping_scores$aupr[[newmapping_id]]
    
    print(cellid)
    print(newauc)
    print(newaupr)
    print("---")
    
    mapping = newmapping_scores$newmapping[[newmapping_id]]
    
    stepid = stepid + 1
    
    steps[[length(steps) + 1]] = tibble(original=list(1:length(mapping)), mapping=list(mapping), stepid=stepid, auc=newauc, aupr=newaupr)
  }
}
steps = bind_rows(steps) %>% mutate(stepid=1:length(steps))
p = steps %>% unnest() %>% ggplot(aes(original, mapping, frame=stepid)) + geom_point()
#gg_animate(p, interval=0.1)
Eregulators2 = Emrna[,order(phenoData(Emrna)$time)[mapping]]

plotdata = bind_rows(
  tibble(original=1:length(mapping), mapping = mapping, mapper="cellshift"), 
  tibble(original=rank(phenoData(Emrna)$time), mapping=rank(phenoData(Emrna)$time)[shift_expression(Emrna, delay, T)], mapper="timeshift"),
  tibble(original=1:length(mapping), mapping=1:length(mapping), mapper="original")
)
plotdata %>% ggplot() + geom_point(aes(original, mapping, color=mapper))

# method 3: http://stats.stackexchange.com/questions/146950/how-can-i-estimate-the-delay-between-two-non-periodic-time-series
# estimate the delay between a known regulator and target, do this for every regulator and target and ...?
# no!!! this assumes the delay is constant across the time series
# maybe do some dynamic time warping for every regulator->target

cellid = 50
mapping[[cellid+1]]
mapping[[cellid]]

steps[[1]]

steps[[length(steps)]]$mapping %>% .[[1]] %>%  plot

bind_rows()

bind_rows(steps)

stepdata = lapply(1:length(steps), function(i) tibble(original=1:length(mapping), mapping=steps[[i]], stepid=i)) %>% bind_rows
library(gganimate)
stepdata %>% ggplot() + geom_line(aes(original, mapping))
