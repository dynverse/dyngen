library(GillespieSSA)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(pheatmap)
library(parallel)
library(PRISM)

G = c(1:100)
tfs = G[1:20]
target2tfs <- lapply(G, function(g) {
  if(g%in%tfs) {
    return(double())
  } else {
    return(sample(tfs, 1+rbinom(1, length(tfs)-1, 0.05)))
  }
})
names(target2tfs) = G

#source("bin/ba_network.R")
#G = c(1:50)
#ba.network <- generate.ba(amnt.nodes = length(G), amnt.edges = 100)
#ba.net.df <- network.to.df(ba.network) # indien nodig
#tfs = G[G %in% ba.net.df$i]
#target2tfs = setNames(ba.network$neighbours, G)

## generating reaction formulas (the probability that a reaction occurs in [t, t dt])
kterms = c()
production = mapply(function(g, tfs) {
  if (length(tfs) > 0) {
    inputs = sapply(tfs, function(tf) paste0("(x", tf, "/k", tf, "g", g, ")"))
    
    kterms <<- c(kterms, sapply(tfs, function(tf) paste0("k", tf, "g", g)))
    
    up = paste0("(a0g", g, "+", paste("a1", inputs, sep="*", collapse="+"), ")")
    down = paste0("(1+", paste0(inputs, collapse="+"), ")")
  } else {
    up = paste0("(a0g", g, ")")
    down = "1"
  }
  return(paste0("r", g, "*", up, "/", down))
}, tfs=target2tfs, g=names(target2tfs))

decay = mapply(function(g, tfs) {
  return(paste0("d", g, "*x" , g))
}, tfs=target2tfs, g=names(target2tfs))

formulas = c(production, decay)

## generating the (initial) parameters of the system

R = rlnorm(length(G), log(10))
names(R) = paste0("r", G)
D = rep(1,length(G))
names(D) = paste0("d", G)
K = sapply(kterms, function(kterm) R[paste0("r",str_replace(kterm, "k(\\d*)g\\d*", "\\1"))]/2)
names(K) = kterms
A0 = rep(0.1, length(G))
names(A0) = paste0("a0g", G)

nu <- cbind2(diag(nrow=length(G), ncol=length(G)), -diag(nrow=length(G), ncol=length(G)))
X0 = round(R/2, 0)
names(X0) = paste0("x", G)

## simulation

# function to postprocess the ssa output
process_ssa = function(out, starttime=0) {
  tf = out$args$tf
  rownames(out$data) = NULL
  data = as.data.frame(out$data)
  lastrow = tail(data, n=1)
  lastrow$V1 = tf
  data[nrow(data)+1,] <- lastrow
  data <- data[data$V1<=tf,]
  times = data$V1 + starttime
  data = as.matrix(data[,-1])
  
  rownames(data) = NULL
  
  return(list(times=times, expression=data))
}

## Simulate multiple single cells following the same linear trajectory

# we will change the basal expression level of transcription factors during the simulation
# we need to determine these perturbations beforehand as they will be the same for every cell
# initial basal expression levels
nruns = 20
A0[] = 0.05
A0[sample(tfs, 1+rbinom(1, length(tfs)-1, 0.5))] = 1
params = c(a1=1, R, D, K, A0)

A0s = lapply(c(1:(nruns+1)), function(i) {
  if(i == 1) {
    return(A0)
  } else {
    A0[sample(tfs, 2)] <<- 1
    A0[sample(tfs, 2)] <<- 0.05
    return(A0)
  }
})

simulate_cell = function(timeofsampling=NULL) {
  require(fastgssa)
  
  totaltime = 0
  tf = 5
  burntf = 4
  
  A0 = A0s[[1]]
  params = c(params, A0)
  
  out <- ssa(X0,formulas,nu,params,tf=tf, method = "D")
  output = process_ssa(out)
  expression = output$expression
  times = output$times
  burnin = times > burntf
  expression = output$expression[burnin,]
  times = times[burnin] - burntf
  
  # other runs: perturb TF activity
  tf=1
  for(i in c(1:20)) {
    X0 = tail(expression, n=1)[1,]
    
    A0 = A0s[[i+1]]
    params = c(params, A0)
    
    out <- fastgssa::ssa(X0,formulas,nu,params,tf=tf, method="D")
    output <- process_ssa(out, last(times))
    expression<-methods::rbind2(expression, output$expression)
    times <- c(times, output$times)
    
    if(!is.null(timeofsampling) && last(times) > timeofsampling) {
        break
    }
  }
  
  # return either the full expression matrix, or the expression at timeofsampling
  if(is.null(timeofsampling)) {
    rownames(expression) = c(1:nrow(expression))
    return(list(expression=expression, times=times))
  } else {
    return(expression[findInterval(timeofsampling, times),])
  }
}
celltimes = runif(200, 0, 20)

library(profvis)
profvis({simulate_cell(1)})

#cells = mclapply(celltimes, simulate_cell, mc.cores=8)
cells = qsub.lapply(celltimes, simulate_cell)

E = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
E = E[,apply(E, 2, sd) > 0]

pheatmap(t(E), scale="row", cluster_cols=F, cluster_rows=T, clustering_distance_rows = "correlation", annotation_col = data.frame(time=celltimes, row.names = rownames(E)))

library(ggplot2)
datadf = melt(E)
colnames(datadf) = c("time", "gene", "count")

datadf = datadf[datadf$gene %in% c("x1", "x2", "x3", "x4", "x5"),]

ggplot(datadf) + geom_line(aes(time, count, group=gene, color=gene))


library(SCORPIUS)
space = reduce.dimensionality(correlation.distance(E),ndim = 2)
trajectory = infer.trajectory(space)
draw.trajectory.plot(space, celltimes, trajectory$final.path) + scale_colour_distiller(palette = "RdYlBu")
rownames(E) = c(1:nrow(E))
draw.trajectory.heatmap(E, trajectory$time, as.factor(cut(celltimes, breaks=length(celltimes), labels=F)))

draw.trajectory.heatmap(E, celltimes, as.factor(cut(celltimes, breaks=99, labels=F)), show.labels.row = T)

