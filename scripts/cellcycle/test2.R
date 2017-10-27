library(tidyverse)
library(magrittr)
library(pheatmap)

## cell cycle model
B <- function(A1, A2, A3, A4) glue::glue("{A2} - {A1} + ({A3}) * ({A2}) + ({A4}) * ({A1})")
G <- function(A1, A2, A3, A4) glue::glue("2 * ({A4}) * ({A1}) / ({B(A1, A2, A3, A4)} + sqrt(({B(A1, A2, A3, A4)})^2 - 4 * ({A2} - {A1}) * ({A4}) * ({A1})))")
G <- function(A1, A2, A3, A4) glue::glue("2 * ({A4}) * ({A1}) / ({B(A1, A2, A3, A4)})")
ref1 <- c(
  cycd = "cycD0 * mass",
  vatf = "katfp + katfp * actcyca + katfppp * actcyce + katfpppp * cycd",
  tfe = G("vatf", "kitfp+kitfpp * actcycb + kitfppp * actcyca", "jatf", "jitf"),
  vde = "kdep + kdepp * actcyce + kdeppp * actcyca + kdepppp * actcycb",
  vda = "kdap + kdapp * cdc20a",
  cdc14 = "cdc20a",
  tfi = G("kafi * cdc14", "kifip + kifipp * actcycb", "jafi", "jifi"),
  vsi = "ksip + ksipp * tfi",
  vdi = "(kdip + kdipp * actcyca+kdippp * actcycb + kdipppp * actcyce + kdippppp * cycd)/(1+(cdc14/j14di))",
  cdk1cycbcki = "cycb - actcycb - prempf",
  cdk1pcycb = "cycb - actcycb - trib",
  tria = "cyca - actcyca",
  trie = "cyce - actcyce",
  freecki = "cki-trib-tria-trie",
  tfb = G("kafb * actcycb", "kifb", "jafb", "jifb"),
  vsb = "ksbp + ksbpp * tfb",
  vdb = "kdbp + kdbpp * cdh1 + kdbppp * cdc20a",
  wee1 = G("kaweep + kaweepp * cdc14", "kiwee * actcycb", "jawee", "jiwee"),
  vwee = "kweep + kweepp * wee1",
  cdc25 = G("ka25 * actcycb", "ki25p + ki25pp * cdc14", "ja25", "ji25"),
  v25 = "k25p + k25pp * cdc25"
) %>% map_chr(~paste0("(", ., ")"))

# itself
ref1 <- map_chr(ref1, function(x) {
  for (name in names(ref1)) {
    x <- gsub(name, ref1[[name]], x)
  }
  
  x
})
ref1 <- map_chr(ref1, function(x) {
  for (name in names(ref1)) {
    x <- gsub(name, ref1[[name]], x)
  }
  
  x
})

formulae.strings <- c(
  dactcyca = "(ksap + ksapp * tfe) * mass + (vdi + kdia) * tria - (vda + kasa * freecki) * actcyca",
  dactcycb = "vsb * mass + v25 * (cycb - trib - actcycb) + (kdib + vdi) * (cycb - prempf - actcycb) - (vdb + vwee + kasb * freecki) * actcycb",
  dactcyce = "(ksep + ksepp * tfe) * mass + (vdi + kdie) * trie- (vde + kase * freecki) * actcyce",
  dcyca = "(ksap + ksapp * tfe) * mass - vda * cyca",
  dcycb = "vsb * mass - vdb * cycb",
  dcyce = "(ksep + ksepp * tfe) * mass - vde * cyce",
  dcdh1 = "(kah1p + kah1pp * cdc14) * (1-cdh1) / (jah1 + 1 - cdh1) - ((kih1p + kih1pp * actcyca + kih1ppp * actcycb + kih1pppp*actcyce+kih1ppppp * cycd) * cdh1)/(jih1 + cdh1)",
  dcki = "vsi - vdi * cki",
  dtrib = "kasb * (cycb - trib) * freecki - (kdib + vdb + vdi) * trib",
  dprempf = "vwee * (cycb - prempf) - (v25 + vdb) * prempf",
  apcp = "(kaapc * actcycb * (1-apcp))/(jaapc + 1 - apcp) - (kiapc * apcp)/(jiapc + apcp)",
  dcdc20a = "(ka20 * apcp * (cdc20t - cdc20a))/(ja20 + cdc20t - cdc20a) - ki20 / (ji20 + cdc20a) + kd20 * cdc20a",
  dcdc20t = "(ks20p + ks20pp * actcycb^n)/(j20^n + actcycb^n) - kd20 * cdc20t",
  dmass = "mu * mass"
  # dmass = "ifelse(actcycb < 0.3, -mass/2, mu * mass)"
)

formulae.strings <- map_chr(formulae.strings, function(f) {
  for (name in names(ref1)) {
    f <- gsub(name, ref1[[name]], f)
  }
  
  f
})

params <- c(
  mu = log(2)/(10),
  j20 = 100,
  ja20 = 0.005,
  jaapc = 0.01,
  jafb = 0.1,
  jah1=0.01,
  jatf = 0.01,
  ji20=0.005,
  jiapc = 0.01,
  jifb = 0.1,
  jih1=0.01,
  jitf = 0.01,
  ka20=0.0833,
  kaapc = 0.0117,
  kafb = 0.167,
  kah1p = 0.175,
  kah1pp=2.33,
  kasa=16.7,
  kase=16.7,
  katfpp=0.05,
  katfppp=0.0833,
  katfpppp=0.055,
  kd20=0.025,
  kdap=0.000333,
  kdapp=0.333,
  kdbp=0.000833,
  kdbpp=0.333,
  kdbppp=0.0167,
  kdep=0.00167,
  kdepp=0.0167,
  kdeppp=0.167,
  kdepppp=0.167,
  kdia=0.167,
  kdie=0.167,
  kdip=0.167,
  kdipp=0.833,
  kdippp=1.67,
  kdipppp=0.833,
  ki20=0.417,
  kiapc=0.03,
  kifb=0.0167,
  kih1pp=0.2,
  kih1ppp=0.667,
  kitfp=0.0417,
  kitfpp=0.0167,
  kitfppp=0.0167,
  ks20pp=2.5,
  ksap=0,
  ksapp=0.00417,
  ksbp=0.00167,
  ksbpp=0.005,
  ksep=0.00133,
  ksepp=0.05,
  ksip=0.333,
  n=1,
  cycD0=0.5
)

params[c("ja25", "jafi", "jawee", "ji25", "jifi", "jiwee", "j14di", "k25p", "k25pp", "ka25", "kafi", "kasb", "katfp", "kaweep", "kaweepp", "kdappp", "kdib", "kdippppp", "ki25p", "ki25pp", "kifip", "kifipp", "kih1p", "kih1pppp", "kih1ppppp", "kiwee", "ks20p", "ksap", "ksipp", "kweep", "kweepp")] = 0.00001

molecules <- gsub("d(.*)", "\\1", names(formulae.strings))
nus <- diag(length(formulae.strings)) %>% set_colnames(names(formulae.strings)) %>% set_rownames(molecules)

initial.state <- runif(length(molecules)) %>% set_names(molecules)
initial.state <- rep(0, length(molecules)) %>% set_names(molecules)
initial.state[["mass"]] <- 1

system <- list(
  params = params,
  formulae.strings = formulae.strings,
  nus = nus,
  initial.state = initial.state
)

extra_func <- function(x, xprev) {
  if(xprev[["actcycb"]] >= 0.3) {
    if(x[["actcycb"]] < 0.3) {
      x[["mass"]] = x[["mass"]]/2
    }
  }
  x
}

initial.state <- runif(length(molecules)) %>% set_names(molecules)
final.time = 240
out <- fastgssa::ssa(
  system$initial.state, 
  system$formulae.strings,
  system$nus,
  final.time, 
  system$params, 
  method=fastgssa::ssa.em(h=final.time/10000, noise_strength=0), 
  recalculate.all = TRUE, 
  stop.on.negstate =FALSE, 
  stop.on.propensity=FALSE,
  extra.functions=list(extra_func)
)
plot(out$timeseries[, "mass"], type="l")

plot(out$timeseries[, "actcycb"], type="l")
plot(out$timeseries[, 3])
out$timeseries[, molecules] %>% pheatmap(cluster_cols=F, cluster_rows=F, scale="column")


exp <- out$timeseries#[, molecules]
