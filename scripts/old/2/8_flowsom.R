library(FlowSOM)

Efiltered =  Eprot[apply(exprs(Eprot), 1, sd) > 0,apply(exprs(Eprot), 2, sd) > 0]
ff = new("flowFrame", exprs=exprs(Efiltered) %>% t)
result = FlowSOM(ff, colsToUse = featureNames(Efiltered), nClus=4)

PlotStars(result[[1]], view="MST")
