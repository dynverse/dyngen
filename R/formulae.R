AbstractFormula <- setClass(
  "AbstractFormula",
  slots = c(
    type = "character",
    name = "character",
    arguments = "list",
    string = "character"
  ),
  validity = function(object) {
    if (! object@type %in% c("infix.fun", "prefix.fun", "var", "con")) {
      return("type must be one of: infix.fun, prefix.fun, var or con")
    }
    if (object@name == "") return("name must be at least length 1")
    return(TRUE)
  }
)

setMethod("show", "AbstractFormula", function(object) cat("Formula: ", slot(object, "string"), "\n", sep=""))


#' Extract all variables from a formula
#'
#' @param formula An S4 class \code{AbstractFormula}
#'
#' @return A list of S4 class objects
#' @export
#'
#' @examples
#' extract.variables(fvar("c", 1) %f*% (fcon(1) - fvar("c", 2)))
extract.variables <- function(formula) {
  if (formula@type == "var") {
    formula
  } else if (formula@type %in% c("infix.fun", "prefix.fun")) {
    unique(unlist(lapply(formula@arguments, extract.variables), recursive = F))
  } else {
    c()
  }
}

#' Create a new variable
#'
#' @param prefix the prefix of the variable
#' @param ... optional arguments for the variable name
#'
#' @return an S4 class \code{AbstractFormula} object
#' @export
#'
#' @examples
#' total <- fvar("total")
#' gene <- fvar("gene", 1)
fvar <- function(prefix, ...) {
  arguments <- list(...)
  AbstractFormula(type = "var", name = prefix, arguments = arguments, string = paste(c(prefix, arguments), collapse="_"))
}

#' Create a new constant
#'
#' @param value the value of the constant
#'
#' @return an S4 class \code{AbstractFormula} object
#' @export
#'
#' @examples
#' con1 <- fcon(1)
#' con2 <- fcon(2.1)
fcon <- function(value) {
  AbstractFormula(type = "con", name = as.character(value), arguments = list(), string = as.character(value))
}

#' Create a new infix function to be used for the AbstractFormula system
#'
#' @param fun the function as a string
#'
#' @return an infix function which to create S4 class \code{AbstractFormula} objects
#' @export
#'
#' @examples
#' `%f-%` <- finfix("-")
finfix <- function(fun) {
  function(lhs, rhs) {
    if (class(lhs) != "AbstractFormula") lhs <- fcon(lhs)
    if (class(rhs) != "AbstractFormula") rhs <- fcon(rhs)
    AbstractFormula(type = "infix.fun", name = fun, arguments = c(lhs, rhs), string = paste0("(", lhs@string, " ", fun, " ", rhs@string, ")"))
  }
}

#' Create a new prefix function to be used for the AbstractFormula system
#'
#' @param fun the function as a string
#'
#' @return a prefix function which to create S4 class \code{AbstractFormula} objects
#' @export
#'
#' @examples
#' fsum <- fprefix("sum")
fprefix <- function(fun) {
  function(...) {
    arguments <- unlist(list(...))
    for (i in seq_along(arguments)) {
      if (class(arguments[[i]]) != "AbstractFormula") {
        arguments[[i]] <- fcon(arguments[[i]])
      }
    }
    arg.string <- paste(sapply(arguments, function(arg) arg@string), collapse = ", ")
    AbstractFormula(type = "prefix.fun", name = fun, arguments = arguments, string = paste0(fun, "(", arg.string, ")"))
  }
}

#' Divide function
#'
#' @param lhs left hand side S4 class \code{AbstractFormula} object
#' @param rhs right hand side S4 class \code{AbstractFormula} object
#'
#' @return an S4 class \code{AbstractFormula} object
#' @rdname fdiv
#' @export
#'
#' @examples
#' fvar("total") %f/% fcon(10)
`%f/%` <- finfix("/")

#' Multiplication function
#'
#' @param lhs left hand side S4 class \code{AbstractFormula} object
#' @param rhs right hand side S4 class \code{AbstractFormula} object
#'
#' @return an S4 class \code{AbstractFormula} object
#' @rdname dmult
#' @export
#'
#' @examples
#' fvar("total") %f*% fcon(10)
`%f*%` <- finfix("*")

#' Addition function
#'
#' @param lhs left hand side S4 class \code{AbstractFormula} object
#' @param rhs right hand side S4 class \code{AbstractFormula} object
#'
#' @return an S4 class \code{AbstractFormula} object
#' @rdname fadd
#' @export
#'
#' @examples
#' fvar("total") %f+% fcon(10)
`%f+%` <- finfix("+")

#' Subtraction function
#'
#' @param lhs left hand side S4 class \code{AbstractFormula} object
#' @param rhs right hand side S4 class \code{AbstractFormula} object
#'
#' @return an S4 class \code{AbstractFormula} object
#' @rdname fmin
#' @export
#'
#' @examples
#' fvar("total") %f-% fcon(10)
`%f-%` <- finfix("-")

setMethod("+", signature(e1 = "AbstractFormula", e2 = "AbstractFormula"), function(e1, e2) e1 %f+% e2)
setMethod("-", signature(e1 = "AbstractFormula", e2 = "AbstractFormula"), function(e1, e2) e1 %f-% e2)
setMethod("*", signature(e1 = "AbstractFormula", e2 = "AbstractFormula"), function(e1, e2) e1 %f*% e2)
setMethod("/", signature(e1 = "AbstractFormula", e2 = "AbstractFormula"), function(e1, e2) e1 %f/% e2)

setMethod("+", signature(e1 = "AbstractFormula", e2 = "numeric"), function(e1, e2) e1 %f+% e2)
setMethod("-", signature(e1 = "AbstractFormula", e2 = "numeric"), function(e1, e2) e1 %f-% e2)
setMethod("*", signature(e1 = "AbstractFormula", e2 = "numeric"), function(e1, e2) e1 %f*% e2)
setMethod("/", signature(e1 = "AbstractFormula", e2 = "numeric"), function(e1, e2) e1 %f/% e2)

setMethod("+", signature(e1 = "numeric", e2 = "AbstractFormula"), function(e1, e2) e1 %f+% e2)
setMethod("-", signature(e1 = "numeric", e2 = "AbstractFormula"), function(e1, e2) e1 %f-% e2)
setMethod("*", signature(e1 = "numeric", e2 = "AbstractFormula"), function(e1, e2) e1 %f*% e2)
setMethod("/", signature(e1 = "numeric", e2 = "AbstractFormula"), function(e1, e2) e1 %f/% e2)

#' Apply a function over all arguments, then collapse it with another function
#'
#' @param collapse.fun a function which can collapse multiple S4 class \code{AbstractFormula} objects into a new one
#' @param lapply.args arguments for the lapply
#' @param lapply.fun a function for the lapply that will return S4 class \code{AbstractFormula} objects
#'
#' @return an S4 class \code{AbstractFormula} object
#' @export
#'
#' @examples
#' freduce(`%f+%`, c(1,2,3), fcon)
freduce <- function(collapse.fun, lapply.args, lapply.fun) {
  Reduce(collapse.fun, lapply(lapply.args, lapply.fun))
}

#' Sum of AbstractFormula objects
#'
#' @param arguments the S4 class \code{AbstractFormula} objects to be summed
#'
#' @return an S4 class \code{AbstractFormula} object
#' @export
#'
#' @examples
#' fsum(lapply(c(1,2,3), fcon))
#' fsum(fcon(1), fcon(2))
fsum <- fprefix("sum")

#' Product of AbstractFormula objects
#'
#' @param ... the S4 class \code{AbstractFormula} objects to be multiplied
#'
#' @return an S4 class \code{AbstractFormula} object
#' @export
#'
#' @examples
#' fprod(lapply(c(1,2,3), fcon))
#' fprod(fcon(1), fcon(2))
fprod <- fprefix("prod")


#' @title Generate production formulae 
#'
#' @param g the gene in question
#' @param regs the regulators of the gene
#' @param amnt.genes the total number of genes
#'
#' @export
generate.production.formula <- function(g, regs, amnt.genes) {
  if (length(regs) == 0) {
    formula <- fvar("r", g) %f*% fvar("a0g", g)
  } else {
    inputs <- fsum(lapply(regs, function(r) fvar("x", r) %f/% fvar("kg", r, g)))
    numerator <- fvar("a0g", g) %f+% (fvar("a1") %f*% inputs)
    denominator <- fcon(1) %f+% inputs
    formula <- fvar("r", g) %f*% numerator %f/% denominator
  }
  
  nu <- rep(0, amnt.genes)
  nu[[g]] <- 1
  
  list(formula = formula, nu = nu)
}

#' @title Generate decay formulae
#'
#' @param g the gene in question
#' @param amnt.genes the total number of genes
#'
#' @export
generate.decay.formula <- function(g, amnt.genes) {
  formula <- fvar("d", g) %f*% fvar("x", g)
  
  nu <- rep(0, amnt.genes)
  nu[[g]] <- -1
  
  list(formula = formula, nu = nu)
}