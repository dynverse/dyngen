#' @title Generate production formulae 
#'
#' @param g the gene in question
#' @param regs the regulators of the gene
#' @param amnt.genes the total number of genes
#'
#' @export
generate.production.formula <- function(g, regs, amnt.genes) {
  if (length(regs) == 0) {
    variables <- paste0(c("r", "a0g"), g)
    formula <- paste0("r", g, "*a0g", g)
  } else {
    variables <- c(
      paste0(c("r", "a0g", "x"), g),
      paste0("k", regs, "g", g)
    )
    
    inputs <- paste0("(x", regs, "/k", regs, "g", g, ")")
    num <- paste0("(a0g", g, "+", paste("a1", inputs, sep="*", collapse="+"), ")")
    denom <- paste0("(1+", paste0(inputs, collapse="+"), ")")
    formula <- paste0("r", g, "*", num, "/", denom)
  }
  nu <- rep(0, amnt.genes)
  nu[[g]] <- 1
  list(variables = variables, formula = formula, nu = nu)
}

#' @title Generate decay formulae
#'
#' @param g the gene in question
#' @param amnt.genes the total number of genes
#'
#' @export
generate.decay.formula <- function(g, amnt.genes) {
  variables <- paste0(c("d", "x"), g)
  formula <- paste0("d", g, "*x", g)
  nu <- rep(0, amnt.genes)
  nu[[g]] <- -1
  list(variables = variables, formula = formula, nu = nu)
}