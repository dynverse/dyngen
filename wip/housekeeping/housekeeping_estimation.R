#### ADD EXTRA GENES

# Different methods:
# * Estimating from real data
#      Simply a poisson distribution
# * Simulating the network independently
#      Really worth it?
# * Simulating together
#      Really worth it?



add_housekeeping_poisson <- function(expression, geneinfo, ngenes=200, overallaverage = mean(expression)) {
  reference_expression <- (2^SCORPIUS::ginhoux$expression)-1
  meanpoissons <- colMeans(reference_expression) %>% {./mean(reference_expression)*overallaverage}
  
  additional_expression <- purrr::map(sample(meanpoissons, ngenes), ~rpois(nrow(expression), .)) %>% 
    invoke(cbind, .) %>% 
    magrittr::set_colnames(paste0("G", seq_len(ngenes)+ncol(expression)))
  
  geneinfo <- dplyr::bind_rows(geneinfo %>% dplyr::mutate(housekeeping=F), tibble(gene=colnames(additional_expression), housekeeping=T))
  
  list(expression=cbind(expression, additional_expression), geneinfo=geneinfo)
}



library(SCORPIUS)
library(tidyverse)
data(ginhoux)

ginhoux$expression %>% dim


rowMea


subexpr = ginhoux$expression[,apply(ginhoux$expression, 2, sd)<1]

apply(subexpr, 2, mean) %>% hist
apply(subexpr, 2, mean) %>% hist



(apply(ginhoux$expression, 2, sd)/apply(ginhoux$expression, 2, mean)) %>% hist


library(zeroinfl)


library(splatter)

library(tidyverse)

expression <- ginhoux$expression[,apply(ginhoux$expression, 2, sd)<1] %>% as.matrix() %>% reshape2::melt(varnames=c("cell", "gene"), value.name="expression") %>% dplyr::mutate(gene=factor(gene), cell=factor(cell)) %>% dplyr::mutate(expression = round(2^expression)-1)
model = expression %>% dplyr::filter(gene %in% unique(gene)[1:600]) %>% group_by(gene) %>% do(model=zeroinfl(expression~1, .))
