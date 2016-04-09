library(fastgssa)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(pheatmap)
library(parallel)
library(PRISM)

library(dyngen)

totaltime = 40
burntime = 0

formulae.strings = c("u_1*y_1*c", "b_1*k_1*c", "r * (a0+b_1*(a1-a0))", "x_2*d", "x_2*p", "y_2*o")
formulae.strings[[3]] = "r * ifelse(b_1, a1, a0)"
formulae.nus = matrix(
  c(-1, -1, 1,0,0,   
    1, 1, -1,0,0,  
    0,0,0,1,0,   
    0,0,0,-1,0, 
    0,0,0,0,1, 
    0,0,0,0,-1
), nrow = 5, dimnames = list(c("y_1", "u_1", "b_1", "x_2", "y_2"), 1:length(formulae.strings)))

params = c(k_1=10,a0=0.05, a1=1, c=0.1, d=1, r=10, o=0.1, p=0.1) 
# c=crowdiness, see introduction Gronlund 2012 Nat Comm
# you could call it a "burst" parameter
# it resembles how hard it is for a transcription factor to get to the binding site
# we use this parameter to keep the k parameter at concentration of half maximal binding

initial.state = c(y_1=10, u_1=1, b_1=0, x_2=0, y_2=0)

output = simulate_cell()

expression = output$expression
times = output$times

binding = data.frame(expression=rep(expression[,"b_1"], each=2)[-length(times)*2], time=rep(1:length(times), each=2)[-1], var="b_1")
ggplot(binding) + geom_line(aes(time, expression))
pracma::trapz(binding$time, binding$expression)/max(binding$time) # percentage of binding

E = expression[,str_detect(colnames(expression), "x_"),drop=F]
Eprot = expression[,str_detect(colnames(expression), "y_"),drop=F]

plotdata = bind_rows(
  binding,
  melt(E, varnames=c("time", "var"), value.name="expression"),
  melt(Eprot, varnames=c("time", "var"), value.name="expression")
)
plotdata$type = sapply(str_split(plotdata$var, "_"), function(x) x[[1]])
plotdata$gene = sapply(str_split(plotdata$var, "_"), function(x) x[[2]])
plotdata$time = output$times[plotdata$time]
ggplot(plotdata) + geom_line(aes(time, expression, group=var, color=var)) + facet_grid(gene~type, scales = "free")
