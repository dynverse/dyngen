## from http://science.sciencemag.org/content/321/5885/126.full
## see supplementary: http://science.sciencemag.org/content/sci/suppl/2008/07/03/321.5885.126.DC1/Tsai.SOM.pdf

library(tidyverse)
library(magrittr)
library(pheatmap)

formulae.strings <- c(
  dcyclin = "ksynth - kdest * apcact * cyclin - ka * (cdk1tot - cdk1cyclin - cdk1cyclinyp - cdk1cyclinyptp - cdk1cyclintp) * cyclin - kd * cdk1cyclin",
  
  dcdk1cyclin = "ka * (cdk1tot - cdk1cyclin - cdk1cyclinyp - cdk1cyclinyptp - cdk1cyclintp) * cyclin - kd * cdk1cyclin - kdest * apcact * cdk1cyclin - kwee1 * wee1act * cdk1cyclin - kwee1basal * (wee1tot - wee1act) * cdk1cyclin + kcdc25 * cdc25act * cdk1cyclinyp +kcdc25basal * (cdc25tot - cdc25act) * cdk1cyclinyp",
  
  dcdk1cyclinyp = " kwee1 * wee1act * cdk1cyclin + kwee1basal * (wee1tot - wee1act) * cdk1cyclin -
  kcdc25 * cdc25act * cdk1cyclinyp -  kcdc25basal * (cdc25tot - cdc25act) *  cdk1cyclinyp - 
  kcak * cdk1cyclinyp + kpp2c * cdk1cyclinyptp -  kdest * apcact * cdk1cyclinyp",
  
  dcdk1cyclinyptp = "kcak * cdk1cyclinyp -  kpp2c * cdk1cyclinyptp - 
  kcdc25 * cdc25act * cdk1cyclinyptp - kcdc25basal * (cdc25tot - cdc25act) * cdk1cyclinyptp +
  kwee1 * wee1act * cdk1cyclintp + kwee1basal * (wee1tot - wee1act) * cdk1cyclintp - 
  kdest * apcact * cdk1cyclinyptp",
  
  dcdk1cyclintp = "kcdc25 * cdc25act * cdk1cyclinyptp + kcdc25basal * (cdc25tot - cdc25act) * cdk1cyclinyptp - kwee1 * wee1act * cdk1cyclintp -
  kwee1basal * (wee1tot-wee1act) * cdk1cyclintp - kdest * apcact * cdk1cyclintp",
  
  dcdc25act = "kcdc25on * (cdk1cyclintp^ncdc25)/(ec50cdc25^ncdc25 + cdk1cyclintp^ncdc25) * (cdc25tot - cdc25act) - kcdc25off * cdc25act",
  
  dwee1act = "-kwee1off * (cdk1cyclintp^nwee1)/(ec50wee1^nwee1 + cdk1cyclintp^nwee1) * wee1act + kwee1on * (wee1tot - wee1act)",
  
  dplx1act = "kplx1on * (cdk1cyclintp^nplx1)/(ec50plx1^nplx1 + cdk1cyclintp^nplx1) * (plx1tot - plx1act) - kplx1off * plx1act",
  
  dapcact = "kapcon * (plx1act^napc)/(ec50plx1^napc + plx1act^napc) * (apctot - apcact) - kapcoff * apcact"
)

params <- c(
  kdest=0.1,
  ka=0.1,
  kd=0.001,
  kwee1=0.05,
  kcdc25=0.1,
  cdc2tot=230,
  cdc25tot=15,
  wee1tot=15,
  apctot=50,
  plx1tot=50,
  nwee1=4,
  ncdc25=4,
  napc=4,
  nplx1=4,
  ec50plx1=40,
  ec50wee1=40,
  ec50cdc25=40,
  ec50apc=40,
  kcdc25on=1.75,
  kcdc25off=0.2,
  kapcon=1,
  kapcoff=0.15,
  kplx1on=1,
  kplx1off=0.15,
  kwee1on=0.2,
  kwee1off=1.75,
  kcak=0.8,
  kpp2c=0.008,
  ksynth=1,
  cdk1tot=100
)
r = 10
params[["kwee1basal"]] <- params[["kwee1"]]/r
params[["kcdc25basal"]] <- params[["kcdc25"]]/r
params[["cycle_speed"]] <- 100

formulae.strings <- map_chr(formulae.strings, ~glue::glue("cycle_speed * ({.})"))

molecules <- gsub("d(.*)", "\\1", names(formulae.strings))
nus <- diag(length(formulae.strings)) %>% set_colnames(names(formulae.strings)) %>% set_rownames(molecules)

initial.state <- rep(0, length(molecules)) %>% set_names(molecules)
initial.state[["wee1act"]] <- params[["wee1tot"]]

system <- list(
  params = params,
  formulae.strings = formulae.strings,
  nus = nus,
  initial.state = initial.state
)

final.time = 10
out <- fastgssa::ssa(system$initial.state, system$formulae.strings, system$nus, final.time, system$params, method=fastgssa::ssa.em(h=final.time/10000, noise_strength=0), recalculate.all = TRUE, stop.on.negstate =FALSE, stop.on.propensity=FALSE)

plot(out$timeseries[, "cyclin"])
plot(out$timeseries[, "wee1act"])
plot(out$timeseries[, "cdk1cyclin"])
plot(out$timeseries[, "cyclin"], out$timeseries[, "cdk1cyclintp"], type="l")
out$timeseries[, molecules] %>% pheatmap(cluster_cols=F, cluster_rows=F)


source("../dynverse/dynmodular/dimred_wrappers.R")
expression <- out$timeseries[, molecules]
expression <- expression[sort(sample(nrow(expression), 10000)), ]
space <- dimred_ica(expression, 2)

space %>% as.data.frame %>% ggplot() + geom_path(aes(Comp1, Comp2))


expression <- out$timeseries[, molecules]

saveRDS(system, "data/systems/cell_cycle_tsai.rds")
