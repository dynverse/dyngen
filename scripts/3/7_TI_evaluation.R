source("scripts/2/SCORPIUS-paper/code_benchmarks/methods/scorpius.R")
source("scripts/2/SCORPIUS-paper/code_benchmarks/methods/monocle.R")
source("scripts/2/SCORPIUS-paper/code_benchmarks/methods/wanderlust.R")
source("scripts/2/SCORPIUS-paper/code_benchmarks/methods/embeddr.R")
source("scripts/2/SCORPIUS-paper/code_benchmarks/methods/tscan.R")

Efiltered =  Emrna2[apply(exprs(Emrna2), 1, sd) > 0,apply(exprs(Emrna2), 2, sd) > 0]

counts = Efiltered %>% exprs %>% t
group.index = phenoData(Efiltered)$time %>% cut(5, labels=F)

settings = list(
  list(name="scorpius", func=scorpius, params=scorpius.parameters),
  list(name="monocle", func=monocle, params=monocle.parameters),
  list(name="wanderlust", func=wanderlust, params=wanderlust.parameters),
  list(name="embeddr", func=embeddr, params=embeddr.parameters),
  list(name="tscan", func=tscan, params=tscan.parameters)
)

results = lapply(settings, function(setting) {
  params = setting$params
  params$counts = counts
  params$group.index = group.index
  results = do.call(setting$func, params)
  results$name = setting$name
  results
})

scores = tibble(name=map_chr(results, ~.$name), cor=map_dbl(results, ~abs(cor(phenoData(Efiltered)$time, .$value))))

scores %>% ggplot() + geom_bar(aes(name, cor), stat="identity")

result = results[[2]]
