library(tidyverse)
library(magrittr)
library(pheatmap)

## cell cycle model
ref1 <- c(
  v29 ="E2FRB * (K20 * ((actCycD + TriD) * LD + LA * (actCycACDk1 + actCycACdk2) + LB * actCycB + LE * actCycE))",
  v30 ="pE2FRB * (K20 * (LD * (actCycD + TriD) + LA * (actCycACDk1 + actCycACdk2) + LB * actCycB + LE * actCycE))",
  v43 ="RB * (K20 * (LD * (actCycD + TriD) + LA * (actCycACDk1 + actCycACdk2) * LB * actCycB + LE * actCycE))",
  v44 ="ppRB * (K19a * (PP1T - PP1A) + K19 * PP1A)",
  v45 ="K26R * E2FRB",
  v46 ="E2F * (K23a * (actCycACDk1 + actCycACdk2) + K23b * actCycB)",
  v47 ="K22 * pE2F",
  v48 ="K26 * E2F * RB",
  v49 ="K26R * pE2FRB",
  v50 ="K26 * RB * pE2F",
  v51 ="K22 * pE2FRB",
  v52 ="E2FRB * (K23a * (actCycACDk1 + actCycACdk2) + K23b * actCycB)",
  Vatf ="katfpp * (actCycACDk1 + actCycACdk2) + katfppp * actCycE + katfpppp * actCycD",
  Vde ="kdep + kdepp * actCycE + kdeppp * (actCycACDk1 + actCycACdk2) + kdepppp * actCycB",
  Vda ="kdap + kdapp * Cdc20A + kacdh1 * Cdh1",
  TFAB ="(2 * A14 * A11) / ((A12 - A11 + A13 * A12 + A14 * A11) + sqrt((A12 - A11 + A13 * A12 + A14 * A11)^2 - 4 * (A12 - A11) * A13 * A11))",
  Vsi ="ksip",
  Vsb ="ksbp + ksbpp * TFAB + ksbppp * actCycB + ksbppp * E2F",
  Vdb ="kdbp + kdbpp * Cdh1 + kdbppp * Cdc20A",
  Wee1 ="(2 * A24 * A21) / ((A22 - A21 + A23 * A22 + A24 * A21) + sqrt((A22 - A21 + A23 * A22 + A24 * A21)^2 - 4 * (A22 - A21) * A23 * A21))",
  Vwee ="kweep + kweepp * Wee1",
  Cdc25="(2 * A34 * A31) / ((A32 - A31 + A33 * A32 + A34 * A31) + sqrt((A32 - A31 + A33 * A32 + A34 * A31)^2 - 4 * (A32 - A31) * A33 * A31))",
  V25 ="k25p + k25pp * Cdc25",
  Vdi ="(kdip + kdipp * (actCycACDk1 + actCycACdk2) + kdippp * actCycB + kdipppp * actCycE)",
  TriE ="cycE - actCycE",
  freeCKI ="CKI - TriA - TriE - TriD",
  CdkCycBCKI ="cycB - actCycB - preMPF",
  Cdk1PCycB ="cycB - actCycB",
  PP1A ="(PP1T / K21) * (FE * (actCycACDk1 + actCycACdk2 + actCycE) + FB * actCycB + 1)"
) %>% map_chr(~paste0("(", ., ")"))

ref2 <- c(
  A11 = "kafab * (actCycACDk1 + actCycACdk2)",
  A12 = "kifb",
  A13 = "Jafb",
  A14 = "Jifb",
  A21 = "kaweep",
  A22 = "kiwee * (actCycACDk1 + actCycACdk2) + kiweeb * actCycB",
  A23 = "Jawee",
  A24 = "Jiwee",
  A31 = "ka25 * actCycB",
  A32 = "ki25p",
  A33 = "Ja25",
  A34 = "Ji25"
)

ref1 <- map_chr(ref1, function(ref1) {
  for (name in names(ref2)) {
    ref1 <- gsub(name, ref2[[name]], ref1)
  }
  
  ref1
})

# itself
ref1 <- map_chr(ref1, function(x) {
  for (name in names(ref1)) {
    x <- gsub(name, ref1[[name]], x)
  }
  
  x
})

formulae.strings[["dE2F"]]

formulae.strings <- c(
  dERG = "k15/(1+(DRG/J15)^2) - k16 * ERG",
  dDRG = "k17p * ERG + (k17 * (DRG/J17)^2)/(1+(DRG/J17)^2) - K18  * DRG",
  dppRB = "v29 + v30 + v43 - v44",
  dE2F = "v29 + v45 + v47 - v46 - v48 + ke2f * E2F * mass - kde2fcdc20 * E2F * Cdc20A - kde2fcdh1 * E2F * Cdh1",
  dpE2F = "v30 + v49 + v46 - v47 - v50 - kde2fcdc20 * pE2F * Cdc20A - kde2fcdh1 * pE2F * Cdh1",
  dRB = "v44 + v45 + v49 - v48 - v50 - v43",
  dE2FRB = "v51 + v48 - v52 - v29 - v45",
  dpE2FRB = "v52 + v50 - v51 - v30 - v49",
  dactCycD = "k9 * DRG + Vdi * TriD + k24r * TriD - k24 * actCycD * freeCKI - k10 * actCycD",
  dTriD = "k24 * actCycD * freeCKI - k24r * TriD - Vdi * TriD - k10 * TriD",
  dactCycACDk1 = "a1frac * ((ksap + ksapp * E2F + ksappp * TFAB) * mass * 2 + (Vdi + kdia) * TriA) - (Vda + kasa * freeCKI) * actCycACDk1",
  dactCycACdk2 = "(1 - a1frac) * ((ksap + ksapp * E2F + ksappp * TFAB) * mass * 2 + (Vdi + kdia) * TriA) - (Vda + kasa * freeCKI) * actCycACdk2",
  dactCycB = "Vsb * mass * 2 + V25 * (cycB - actCycB) - (Vdb + Vwee) * actCycB",
  dactCycE = "(ksep + ksepp * E2F) * mass * 2 + (Vdi + kdie) * TriE - (Vde + kase * freeCKI) * actCycE",
  dcycA = "(ksap + ksapp * E2F + ksappp * TFAB) * mass * 2 - Vda * cycA",
  dcycB = "Vsb * mass * 2 - Vdb * cycB",
  dcycE = "(ksep + ksepp * E2F) * mass * 2 - Vde * cycE",
  dCKI = "Vsi - Vdi * CKI",
  dCdh1 = "((kah1p + kah1pp * Cdc20A) * (1-Cdh1))/(Jah1 + 1 - Cdh1) - ((kih1pp * (actCycACDk1 + actCycACdk2) + kih1ppp * actCycB) * Cdh1)/(Jih1 + Cdh1)",
  dpreMPF = "Vwee * (cycB - preMPF) - (V25 + Vdb) * preMPF",
  dTriA = "kasa * (cycA - TriA) * freeCKI - (kdia + Vda + Vdi) * TriA",
  dAPCP = "(kaAPC * actCycB * (1-APCP)) / (JaAPC + 1 - APCP) - (kiAPC * APCP) / (JiAPC + APCP)",
  dCdc20A = "(ka20 * APCP * (Cdc20T - Cdc20A)) / (Ja20 + Cdc20T - Cdc20A) - (ki20 / (Ji20 + Cdc20A) + kd20) * Cdc20A",
  dCdc20T = "(ks20pp * actCycB) / (J20 + actCycB) - kd20 * Cdc20T",
  dmass = "u * mass"
)


formulae.strings <- map_chr(formulae.strings, function(f) {
  for (name in names(ref1)) {
    f <- gsub(name, ref1[[name]], f)
  }
  
  f
})

params <- c(
  a1frac=0.081283,
  FB=2,
  FE=25,
  J15=0.1,
  J17=0.3,
  J20=100,
  Ja20=0.005,
  Ja25=0.005,
  JaAPC=0.01,
  Jafb = 0.01,
  Jah1=0.15,
  Jatf=0.01,
  Jawee=0.05,
  Jaweeb=0.05,
  Ji20=0.005,
  Ji25=0.031623,
  JiAPC=0.001,
  Jifb=0.001,
  Jih1=0.01,
  Jitf=0.01,
  Jiwee=0.05,
  k10=88.175,
  k15=5.2905,
  k16=44.0875,
  k17=2645.25,
  k17p=2.64525,
  K18=176.35,
  K19=35.27,
  K19a=440.875,
  K20=176.35,
  K21=1,
  K22=3.527,
  K23a=0.17635,
  K23b=1.7635,
  k24=1763.5,
  k24r=176.35,
  k25p=61.474,
  k25pp=30515.96,
  K26=17635,
  K26R=35.27,
  k9=45.851,
  ka20=292.669,
  ka25=8.85277,
  kaAPC=2.33401,
  kacdh1=264.525,
  kafab=0.296268,
  kah1p=155.8708,
  kah1pp=176350,
  kasa = 19733.57,
  kase=19733.57,
  katfpp=58.70692,
  katfppp=97.80724,
  katfpppp=77.63932,
  kaweep=13.8188,
  kd20=17.635,
  kdap=.516094,
  kdapp=2645.25,
  kdbp=0.853181,
  kdbpp=176.35,
  kdbppp=387.97,
  kde2fcdc20=881.75,
  kde2fcdh1=1.7635,
  kdep=1.961012,
  kdepp=1.973357,
  kdeppp=176.35,
  kdepppp=3527,
  kdia=196.0783,
  kdie=196.0783,
  kdip=196.0783,
  kdipp=978.0688,
  kdippp=1960.837,
  kdipppp=978.0688,
  ke2f=4.2324,
  ki20=17.635,
  ki25p=35.27,
  kiAPC=3.862259,
  kifb=9.827456,
  kih1pp=17635,
  kih1ppp=1763.5,
  kitfp=48.96181,
  kitfpp=19.60836,
  kitfppp=19.60836,
  kiwee=0.145,
  kiweeb=5,
  ks20pp=105.81,
  ksap=16.75325,
  ksapp=0.10581,
  ksappp=20.28025,
  ksbp=6.7013,
  ksbpp=15.8715,
  ksbppp=1.7635,
  ksbpppp=0.617225,
  ksep=1.562461,
  ksepp=8.8175,
  ksip=390.9926,
  kweep=234.8312,
  kweepp=17635,
  LA=30,
  LB=0.5,
  LD=3.3,
  LE=10,
  PP1T=1,
  u=0.693937
)

molecules <- gsub("d(.*)", "\\1", names(formulae.strings))
nus <- diag(length(formulae.strings)) %>% set_colnames(names(formulae.strings)) %>% set_rownames(molecules)

initial.state <- runif(length(molecules)) %>% set_names(molecules)
initial.state <- rep(0.001, length(molecules)) %>% set_names(molecules)
initial.state[["actCycD"]] <- 0.1

system <- list(
  params = params,
  formulae.strings = formulae.strings,
  nus = nus,
  initial.state = initial.state
)

final.time = 1
out <- fastgssa::ssa(system$initial.state, system$formulae.strings, system$nus, final.time, system$params, method=fastgssa::ssa.em(h=final.time/1000, noise_strength=0), recalculate.all = TRUE, stop.on.negstate =FALSE, stop.on.propensity=FALSE)

plot(out$timeseries[, 3])
out$timeseries %>% colnames
out$timeseries[, molecules] %>% pheatmap(cluster_cols=F, cluster_rows=F)


exp <- out$timeseries#[, molecules]
