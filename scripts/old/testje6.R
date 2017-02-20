## More complex networks, based on a given network of lineage-determining TFs following a certain trajectory, and which then influence the expression of other TFs

library(fastgssa)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(pheatmap)
library(parallel)
library(PRISM)
qsub.conf = qsub.configuration(exec.before = c("module unload python", "module load python"))
library(readr)
library(Biobase)
library(igraph)

library(dyngen)

source("bin/testje4_functions.R")

#net = read_tsv("data/networks/linear.tsv")
net = read_tsv("data/networks/cycle.tsv")
ldtfs = sort(unique(c(net$from, net$to)))

# source("bin/ba_network.R")
# ba.network <- generate.ba(amnt.nodes = 50, amnt.edges = 50, offset.exponent = 0)
# ba.net.df <- network.to.df(ba.network) # indien nodig
# ba.net.df$effect = 1
# ba.net.df = ba.net.df[!(ba.net.df$to %in% ldtfs),]
# net = bind_rows(net, ba.net.df)
# net = unique(net)

### ONLY RUN IF YOU WANT EXTRA TARGET GENES
subtfs = c(last(ldtfs)+1:10)
net = bind_rows(net, data.frame(from=sample(ldtfs, length(subtfs), replace = T), to=subtfs, effect=1, strength=1, cooperativity=1))

tfs = c(ldtfs, subtfs)

targets = c(last(subtfs)+1:50)
net = bind_rows(net, data.frame(from=sample(tfs, length(targets), replace = T), to=targets, effect=1, strength=1, cooperativity=1))
###

G = sort(unique(union(net$from, net$to)))
tfs = sort(unique(net$from))

graph = graph_from_data_frame(net)
layout <- layout_with_fr(graph)
plot.igraph(graph, edge.color = c("red", "blue")[as.numeric(factor(net$effect, levels = c(-1,1)))], layout=layout, vertex.size = c(5,20)[as.numeric(factor(G %in% ldtfs, levels = c(F, T)))], edge.arrow.size=0.5, vertex.label=lapply(G, function(g) ifelse(g%in%ldtfs, g, "")), edge.loop.angle=0.1)

barplot(degree(graph, mode="in"))

pheatmap(sapply(G, function(i) sapply(G, function(j) .jacc(net$from[net$to==i], net$from[net$to==j]))))

## generating reaction formulas (the probability that a reaction occurs in [t, t dt])
formulae = list()
nus.changes = list()
add.formula = function(formula, nu) {formulae <<- c(formulae, formula);nus.changes <<- c(nus.changes, list(nu))}
variables = list() # store every variable including some extra info for later use
vargroups = list() # store every variable in groups, for example if we want all k values later on
add.variable = function(variable, ...) {
  varname = variable@string
  info = list(name=varname, ...)
  variables[[varname]] <<- info
  
  group = str_replace(varname, "^(.*?)[_$].*", "\\1")
  if (group %in% c("b", "u", "x", "y", "k", "r", "d", "p", "q", "a0", "a")) {
    vargroups[[group]] <<- c(vargroups[[group]], varname)
  }
  return(variable)
}
get.nu = function(higher=NULL, lower=NULL) {setNames(c(rep(1, length(higher)), rep(-1, length(lower))), lapply(c(higher, lower), function(x) {x@string}))}

kg = add.variable(fvar("kg")) # global tf binding speed (influences burstiness of transcriptional regulation)
rg = add.variable(fvar("rg")) # global transcription rate
dg = add.variable(fvar("dg")) # global RNA degradation rate
pg = add.variable(fvar("pg")) # global translation rate
qg = add.variable(fvar("qg")) # global protein degradation rate
for (target in G) {
  x = add.variable(fvar("x", target), gene=target)
  y = add.variable(fvar("y", target), gene=target)
  r = add.variable(fvar("r", target), gene=target)
  d = add.variable(fvar("d", target), gene=target)
  q = add.variable(fvar("q", target), gene=target)
  p = add.variable(fvar("p", target), gene=target)
  a0 = add.variable(fvar("a0", target), gene=target)
  
  subnet = net[net$to == target,]
  regs = subnet$from
  effects = subnet$effect
  
  boundvars = list()
  
  # regulator binding sites
  if (length(regs) > 0) {
    for (i in seq_len(length(regs))) {
      regulator = regs[[i]]
      effect = effects[[i]]
      strength = subnet[i,]$strength
      cooperativity = subnet[i,]$cooperativity
      
      u = add.variable(fvar("u", target, regulator), target=target, regulator=regulator)
      b = add.variable(fvar("b", target, regulator), target=target, regulator=regulator)
      k = add.variable(fvar("k", target, regulator), target=target, regulator=regulator, strength=strength)
      
      boundvars = c(boundvars, b)
      
      y_regulator = fvar("y", regulator)
      
      # binding
      if (cooperativity == 2) {
        formula = u * y_regulator * y_regulator * kg
      }else {
        formula = u * y_regulator * kg
      }
      
      nu = get.nu(b, list(y_regulator, u))
      add.formula(formula, nu)
      
      # deassociation
      formula = b * k * kg
      nu = get.nu(list(y_regulator, u), b)
      add.formula(formula, nu)
      
      # add to decision tree
      a = add.variable(fvar("a", target, regulator), target=target, regulator=regulator, effect=effect)
      
      # quick and dirty decision tree, the last interaction will have priority
      if(i == 1) {
        decisiontree = fcon(paste0("ifelse(", b@string, ",", a@string, ",", a0@string, ")"))
      } else {
        decisiontree = fcon(paste0("ifelse(", b@string, ",", a@string, ",", decisiontree@string, ")"))
      }
    }
  } else {
    decisiontree = 1
  }
  
  # mRNA production
  
  formula = rg * r * decisiontree
  nu = get.nu(x)
  add.formula(formula, nu)
  
  # mRNA degradation
  formula = dg * d * x
  nu = get.nu(lower=x)
  add.formula(formula, nu)
  
  if(target %in% net$from) {
    # protein production
    formula = pg * p * x
    nu = get.nu(y)
    add.formula(formula, nu)
    
    # protein degradation
    formula = qg * q * y
    nu = get.nu(lower=y)
    add.formula(formula, nu)
  }
}

formulae.strings <- sapply(formulae, function(fl) fl@string)

## generating the (initial) parameters of the system

R = sapply(vargroups$r, function(r) 20)
D = sapply(vargroups$d, function(d) 1)

K = sapply(vargroups$k, function(k) unname((R/D)[paste0("r_",str_replace(k, "k_\\d*_(\\d*)", "\\1"))]/2/variables[[k]]$strength))

P = sapply(vargroups$p, function(p) 1)
Q = sapply(vargroups$q, function(q) 1)

A0 = sapply(vargroups$a0, function(a0) 0)
A0[[1]] = 1
A0[[4]] = 1
#A0[paste0("a0_", ldtfs)] = 1 # IF CYCLE

A = sapply(vargroups$a, function(a) ifelse(variables[[a]]$effect==1, 1, 0))

initial.state = numeric()
initial.state[vargroups$x] = 0
initial.state[vargroups$y] = 0
initial.state[vargroups$u] = 1
initial.state[vargroups$b] = 0
molecules = names(initial.state)

params = c(R, D, K, P, Q, A0,A, rg=1,dg=1,kg=1,pg=10,qg=10, a1=1)

formulae.nus = sapply(nus.changes, function(nu.changes) {
  nu = sapply(molecules, function(r) 0)
  nu[names(nu.changes)] = nu.changes
  nu
})

## Simulate multiple single cells following the same trajectory

totaltime = 50
burntime = 0

celltimes = runif(400, 0, totaltime)

simulate_cell()

cells = mclapply(celltimes, simulate_cell, mc.cores=8)

cells = qsub.lapply(celltimes, simulate_cell, qsub.config = qsub.conf)

expression = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
E = ExpressionSet(t(expression[,str_detect(colnames(expression), "x_")]), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
#E = apply(E, 2, function(x) {sapply(x, function(y) {rbinom(1, y, 0.3)})})
#E = E[,apply(E, 2, sd) > 0]

Eprot = ExpressionSet(t(expression[,str_detect(colnames(expression), "y_")]), AnnotatedDataFrame(data.frame(time=celltimes, row.names = 1:length(celltimes))))
#Eprot = Eprot[,apply(Eprot, 2, sd) > 0]

pheatmap(exprs(E)[,order(phenoData(E)$time)], scale="row", cluster_cols=F, cluster_rows=T, clustering_distance_rows = "correlation", annotation_col = data.frame(time=phenoData(E)$time, row.names = sampleNames(E)))
pheatmap(exprs(Eprot)[,order(phenoData(Eprot)$time)], scale="row", cluster_cols=F, cluster_rows=F, clustering_distance_rows = "correlation", annotation_col = data.frame(time=phenoData(Eprot)$time, row.names = sampleNames(Eprot)))

plot(exprs(E["x_6",]), exprs(E["x_5",]))

library(ggplot2)
datadf = melt(Eprot[order(celltimes),])
colnames(datadf) = c("time", "gene", "count")

ggplot(datadf) + geom_line(aes(time, count, group=gene, color=gene))

# trajectory
library(SCORPIUS)
Efiltered =  E[apply(exprs(E), 1, sd) > 0,apply(exprs(E), 2, sd) > 0]
space = reduce.dimensionality(correlation.distance(t(exprs(Efiltered))),ndim = 3)

space = tsne::tsne(as.dist(correlation.distance(t(exprs(Efiltered)))), k = 3)
colnames(space) = paste0("Comp", 1:ncol(space))

space = diffusionMap::diffuse(as.dist(correlation.distance(t(exprs((Efiltered))))), neigen=2)
space = space$X[,c(1:2)]
colnames(space) = paste0("Comp", 1:ncol(space))

rgl::plot3d(space, col=phenoData(Efiltered)$time)

plotdata = as.data.frame(space)
plotdata$progression = phenoData(Efiltered)$time
ggplot(plotdata) + geom_point(aes(Comp1, Comp2, color=progression)) #+  scale_colour_distiller(palette = "RdYlBu") + theme_classic()

trajectory = infer.trajectory(space)
draw.trajectory.plot(space, phenoData(Efiltered)$time, trajectory$final.path) + scale_colour_distiller(palette = "RdYlBu") + theme_classic()
rownames(E) = c(1:nrow(E))
draw.trajectory.heatmap(t(exprs(Efiltered)), trajectory$time, as.factor(cut(phenoData(Efiltered)$time, breaks=length(phenoData(Efiltered)$time), labels=F)), show.labels.row = T)

draw.trajectory.heatmap(t(exprs(Efiltered)), phenoData(Efiltered)$time, as.factor(cut(phenoData(Efiltered)$time, breaks=length(phenoData(Efiltered)$time), labels=F)), show.labels.row = T)

# check a random tf and his targets, for example to check the time lag
tf = tfs[[6]]
print(tf)
targets = net[net$from == tf,]$to
length(targets)
pheatmap(t(expression[order(celltimes),c(paste0("y_", tf), paste0("x_", targets))]), scale="row", cluster_rows=F, cluster_cols = F)

# check protein and mRNA
i = 1
x = paste0("x_",i)
y = paste0("y_",i)
plotdata = bind_rows(
  data.frame(expression=E[order(celltimes), x], g=x, time=sort(celltimes)),
  data.frame(expression=Eprot[order(celltimes), y], g=y, time=sort(celltimes))
)

ggplot(plotdata) + geom_point(aes(time, expression, group=g, color=g))

#

nknots = 8
Esmooth = cbind2(apply(E, 2, function(x) {smooth.spline(celltimes, x, nknots=nknots)$y}), apply(Eprot, 2, function(x) {smooth.spline(celltimes,x, nknots=nknots)$y}))
Esmooth = cbind2(apply(E[order(celltimes), ], 2, smooth, kind="3RS3R"), apply(Eprot[order(celltimes), ], 2, smooth, kind="3RS3R"))
plotdata = melt(Esmooth, value.name="expression")
plotdata = melt(cbind2(E, Eprot), value.name="expression")
colnames(plotdata) = c("time", "var", "expression")
plotdata$time = celltimes[plotdata$time]
plotdata$type = sapply(str_split(plotdata$var, "_"), function(x) x[[1]])
plotdata$gene = sapply(str_split(plotdata$var, "_"), function(x) x[[2]])

ggplot(plotdata, aes(x=time, y=expression)) + 
  geom_point(aes(color=gene, shape=type), alpha=.1)+
  geom_smooth(aes(group=var, color=gene, linetype=type)) + 
  facet_wrap(~gene)

# check multiple individual cells
plotdata = bind_rows(lapply(seq_len(10), function(i) {
  out = simulate_cell()
  data.frame(expression=out$expression[,x], simulation=i, time=out$times)
}))
ggplot(plotdata) + geom_line(aes(time, expression, group=simulation, color=simulation))

##
# one individual simulation
totaltime=50
output = simulate_cell()

filtervars = colnames(output$expression)[sapply(strsplit(colnames(output$expression), "_"), function(x) {x[[2]] %in% ldtfs})]

expression = output$expression[,filtervars]
times = output$times

# % of the time bound
binding = expression[,str_detect(colnames(expression), "b_"),drop=F]
#binding = bind_rows(lapply(colnames(binding), function(g) data.frame(expression=rep(binding[,g], each=2)[-length(times)*2], time=rep(1:length(times), each=2)[-1], var=g))) # step plot
dlply(binding, "var", function(binding) pracma::trapz(binding$time, binding$expression)/max(binding$time)) # percentage of binding

ma <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=1)}
binding = apply(binding, 2, ma, n=50)
binding = melt(binding, varnames = c("time", "var"), value.name = "expression")

ggplot(binding) + geom_line(aes(time, expression, color=var)) + facet_wrap(~var)


E = expression[,str_detect(colnames(expression), "x_"),drop=F]
Eprot = expression[,str_detect(colnames(expression), "y_"),drop=F]
Ebind= expression[,str_detect(colnames(expression), "b_"),drop=F]

plotdata = bind_rows(
  #melt(Ebind, varnames=c("time", "var"), value.name="expression"),
  melt(apply(Ebind, 2, ma, n=1000), varnames = c("time", "var"), value.name = "expression"),
  melt(E, varnames=c("time", "var"), value.name="expression"),
  melt(Eprot, varnames=c("time", "var"), value.name="expression")
)
plotdata$type = sapply(str_split(plotdata$var, "_"), function(x) x[[1]])
plotdata$gene = sapply(str_split(plotdata$var, "_"), function(x) x[[2]])
plotdata$goi = sapply(str_split(plotdata$var, "_"), function(x) ifelse(x[[1]] == "b", x[[3]], x[[2]]))
plotdata$time = output$times[plotdata$time]
ggplot(plotdata) + geom_step(aes(time, expression, group=var, color=goi)) + facet_grid(gene~type, scales = "free_y")
 
