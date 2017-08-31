#' Compute weighted adjacency matrix of inferred network
#'
#' @param data A matrix of observations of the different features. The rows must contain the observations.
#' @param K The choice of number of input genes randomly selected as candidates at each node. Must be \code{"all"} for all input features, \code{"sqrt"} for the square root of all input features (default), or an integer.
#' @param nb.trees The number of trees in ensemble for each target gene (default 1000).
#' @param regulators A set of indices or column names of entities whose observed values regulate the observed values of the targets.
#' @param targets A set of indices or column names of entities whose observed values are regulated by the regulators.
#' @param importance.measure Type of variable importance measure. Must be either \code{"IncNodePurity"} for the importance measure based on decrease of residual sum of squares, or \code{"\%IncMSE"} for importance measure optained by permutation of OOB data.
#' @param seed A random number generator seed for replication of analyses. NULL means the seed is not set.
#' @param trace Output additional information.
#' @param mc.cores The number of threads to use for parallel execution.\code{\link{genie3}}
#'
#' @references Huynh-Thu, V. A. et al. (2010) Inferring Regulatory Networks from Expression Data Using Tree-Based Methods. PLoS ONE.
#'
#' @return The weighted adjacency matrix of inferred network.
#' @export
#'
#' @importFrom randomForest randomForest
#' @importFrom parallel mclapply
#'
#' @examples
#' library(GENIE3)
#' data <- matrix(runif(100*100), ncol=100)
#' true.matrix <- matrix(sample(0:1, 20*100, replace=TRUE, prob=c(.9, .1)), ncol=10)
#' diag(true.matrix) <- 0
#' weights <- genie3(data, regulators=1:20, targets=1:100, mc.cores=8)
#' ranking <- get.ranking(weights)
#' evaluation <- evaluate.ranking(ranking, true.matrix=true.matrix)
#' evaluation$au.score
#'
#' # draw a ROC curve
#' library(ggplot2)
#' ggplot(evaluation$metrics, aes(fpr, rec)) + geom_line() + coord_cartesian(xlim = c(0, 1), ylim=c(0, 1)) + theme_minimal()
#'
#' # draw a PR curve
#' ggplot(evaluation$metrics, aes(rec, prec)) + geom_line() + coord_cartesian(xlim = c(0, 1), ylim=c(0, .1)) + theme_minimal()
#'
#' # Evaluate multiple rankings at the same time
#' weights.cor <- as.matrix(abs(cor(weights[,1:20], weights[,1:100])))
#' ranking.cor <- get.ranking(weights.cor)
#' rankings <- list(GENIE3=ranking, Correlation=ranking.cor)
#' evaluations <- evaluate.multiple.rankings(rankings, true.matrix=true.matrix)
#' evaluations$au.score
#'
#' # draw a ROC curve
#' ggplot(evaluations$metrics, aes(fpr, rec, colour=ranking.name)) + geom_line() + coord_cartesian(xlim = c(0, 1), ylim=c(0, 1)) + theme_minimal()
#'
#' # draw a PR curve
#' ggplot(evaluations$metrics, aes(rec, prec, colour=ranking.name)) + geom_line() + coord_cartesian(xlim = c(0, 1), ylim=c(0, .1)) + theme_minimal()
genie3 <- function(data, regulators=seq_len(ncol(data)), targets=seq_len(ncol(data)), data.regulators=NULL,
                   K="sqrt", nb.trees=1000, importance.measure="IncNodePurity",
                   seed=NULL, trace=TRUE, mc.cores=1) {
  requireNamespace("randomForest")
  
  if (!(is.matrix(data) || is.data.frame(data)) || any(!apply(data, 2, is.finite))) {
    stop("Parameter \"data\" must be a matrix or a data frame consisting of finite values.")
  }
  
  num.samples <- nrow(data)
  num.genes <- ncol(data)
  feature.names <- colnames(data)
  
  # check nb.trees parameter
  if (nb.trees < 1) {
    stop("Parameter \"nb.trees\" must be larger than or equal to 0")
  }
  
  # check regulators parameter
  if (is.null(regulators)) {
    regulators <- seq_len(num.genes)
  } else if (class(regulators) == "character") {
    regulators <- match(regulators, feature.names)
  }
  if ((class(regulators) != "numeric" && class(regulators) != "integer") ||
      !all(regulators %in% seq_len(num.genes)) ||
      length(unique(regulators)) != length(regulators)) {
    stop("Parameter \"regulators\" must either be NULL, or be a subset of either colnames(data) or seq_len(ncol(data)), and may not contain any repeated elements.")
  }
  num.regulators <- length(regulators)
  regulator.names <- feature.names[regulators]
  
  # check targets parameter
  if (is.null(targets)) {
    targets <- seq_len(num.genes)
  } else if (class(targets) == "character") {
    targets <- match(targets, feature.names)
  }
  if ((class(targets) != "numeric" && class(targets) != "integer") ||
      !all(targets %in% seq_len(num.genes)) ||
      length(unique(targets)) != length(targets)) {
    stop("Parameter \"targets\" must either be NULL, or be a subset of either colnames(data) or seq_len(ncol(data)), and may not contain any repeated elements.")
  }
  num.targets <- length(targets)
  target.names <- feature.names[targets]
  
  # check K parameter
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.regulators))
  } else if (K == "all") {
    mtry <- num.regulators-1
  } else {
    stop("Parameter \"K\" must be \"sqrt\", or \"all\", or an integer")
  }
  
  # check importance.measure parameter
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter \"importance.measure\" must be \"IncNodePurity\" or \"%IncMSE\"")
  }
  
  # check seed parameter, use the current seed if none is given.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # check mc.cores parameter
  if (mc.cores == "qsub") {
    requireNamespace("PRISM")
    lapplyfun <- function(X, FUN, ...) PRISM::qsub_lapply(X = X, FUN = FUN, qsub.environment = c("data", "feature.names", "regulators", "nb.trees", "importance.measure", "regulator.names", "num.regulators"))
  } else if (mc.cores < 1) {
    stop("Parameter \"mc.cores\" must be larger than or equal to 0")
  } else if (mc.cores > 1) {
    requireNamespace("parallel")
    lapplyfun <- function(...) parallel::mclapply(..., mc.cores=mc.cores)
  } else {
    lapplyfun <- lapply
  }
  
  if (trace) {
    cat("GENIE3 parameter checks are OK!\n")
    cat("Starting GENIE3 computations.\n")
    flush.console()
  }
  
  # scale the data matrix
  data <- scale(data)
  if (!is.null(data.regulators)) data.regulators = scale(data.regulators)
  
  # compute importances for every target
  output <- do.call("cbind", lapplyfun(targets, function(target.index) {
    if (trace) {
      cat("Computing for target ", ifelse(is.null(feature.names), target.index, feature.names[[target.index]]), "\n", sep="")
      flush.console()
    }
    if (!is.null(data.regulators)) {
      x <- data.regulators[,setdiff(regulators, target.index),drop=F]
    } else {
      x <- data[,setdiff(regulators, target.index),drop=F]
    }
    
    y <- data[,target.index]
    rf <- randomForest::randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=TRUE)
    im <- rf$importance[,importance.measure]
    
    out <- setNames(rep(0, num.regulators), regulator.names)
    out[target.index != regulators] <- im
    out
  }))
  
  rownames(output) <- regulator.names
  colnames(output) <- target.names
  
  output / num.samples
}

#' Make an ordered data frame of regulatory links
#'
#' @param weights A weighted adjacency matrix as returned by \code{\link{genie3}}
#' @param max.links The maximum number of links to output.
#'
#' @return A data frame of links.
#' @import reshape2 dplyr
#' @export
#'
#' @seealso \code{\link{genie3}}
get.ranking <- function(weights, max.links=100000) {
  requireNamespace("reshape2")
  requireNamespace("dplyr")
  
  # This following line is added just to appease R check. When looking at this code, pretend this line doesn't exist.
  regulator <- target <- value <- NULL
  
  ranking <- reshape2::melt(weights, varnames=c("regulator", "target"))                # transform into an expanded data frame
  ranking <- dplyr::filter(ranking, as.character(regulator) != as.character(target))   # remove self loops
  ranking <- dplyr::arrange(ranking, dplyr::desc(value))                               # order by value
  
  if (is.integer(max.links) && nrow(ranking) > max.links) {
    ranking <- dplyr::top_n(ranking, value, max.links)                                 # take top max.links edges
  }
  
  ranking
}

#' Evaluate a network ranking
#'
#' @param ranking the data frame as returned by \code{\link{get.ranking}}. This data frame must contain at least 2 columns, called \code{regulator} and \code{target}.
#' @param true.matrix a matrix with 0's and 1's, representing the golden standard. The rownames and colnames must me the same as the names used in the regulator and target columns in \code{ranking}.
#'
#' @return a list containing 2 items, the ranked evaluation and the area under the curve scores
#' @import ROCR pracma dplyr
#' @export
#'
#' @seealso \code{\link{genie3}}
evaluate.ranking <- function(ranking, true.matrix) {
  requireNamespace("ROCR")
  requireNamespace("pracma")
  requireNamespace("dplyr")
  
  # This following line is added just to appease R check. When looking at this code, pretend this line doesn't exist.
  regulator <- target <- value <- NULL
  
  # get TPs along ranking
  reg.names <- if (class(ranking$regulator) == "integer") seq_len(nrow(true.matrix)) else rownames(true.matrix)
  tar.names <- if (class(ranking$target) == "integer") seq_len(ncol(true.matrix)) else colnames(true.matrix)
  ranking <- dplyr::mutate(ranking,
                           tp=mapply(regulator, target, FUN=function(r, t) if (r %in% reg.names & t %in% tar.names) true.matrix[[r, t]] else 0),
                           value=if("value" %in% colnames(ranking)) value else percent_rank(-seq_len(nrow(ranking)))
  )
  
  evaluate.ranking.direct(ranking$value, ranking$tp)
}

#' Title
#'
#' @param value
#' @param trues
#' @param perf.measures the ROCR performance measures (See \code{\link[ROCR]{performance}}). Must at least contain \code{"fpr"}, \code{"rec"}, \code{"spec"}, and \code{"prec"}.
#'
#' @return
#' @export
#'
#' @examples
evaluate.ranking.direct <- function(value, trues, perf.measures=c("acc", "rec", "prec", "fpr", "spec", "phi", "f")) {
  rocr.pred <- ROCR::prediction(value, trues)
  
  # run different metrics using ROCR
  metrics <- do.call("cbind", lapply(perf.measures, function(p) ROCR::performance(rocr.pred, p, x.measure="cutoff")@y.values[[1]]))
  colnames(metrics) <- perf.measures
  
  # combine into one data frame
  eval.df <- data.frame(
    cutoff = ROCR::performance(rocr.pred, "acc", x.measure="cutoff")@x.values[[1]],
    balanced.acc = rowMeans(metrics[,c("rec", "spec")]),
    metrics
  )
  eval.df[1,c("prec", "phi", "f")] <- 1
  
  # calculate AUs
  auroc <- pracma::trapz(eval.df$fpr, eval.df$rec)
  aupr <- abs(pracma::trapz(eval.df$rec, eval.df$prec))
  #   if (num.gold > 1) {
  #     aupr <- aupr / (1.0 - 1.0 / num.gold)
  #   }
  F1 <- ifelse(auroc + aupr != 0, 2 * auroc * aupr / (auroc + aupr), 0)
  au.df <- data.frame(auroc = auroc, aupr =  aupr, F1 = F1)
  
  # generate output
  list(metrics = eval.df, au.score = au.df)
}

#' Title Evaluate and compare multiple rankings
#'
#' @param rankings a list of rankings as preduced by \code{\link{get.ranking}}. See \code{\link{evaluate.ranking}} for more information.
#' @param true.matrix a matrix with 0's and 1's, representing the golden standard. The rownames and colnames must me the same as the names used in the regulator and target columns in \code{ranking}.
#' @param perf.measures the ROCR performance measures (See \code{\link[ROCR]{performance}}). Must at least contain \code{"fpr"}, \code{"rec"}, \code{"spec"}, and \code{"prec"}.
#'
#' @return a list containing 2 items, the ranked evaluation and the area under the curve scores
#' @export
#'
#' @seealso \code{\link{genie3}}
evaluate.multiple.rankings <- function(rankings, true.matrix, perf.measures=c("acc", "rec", "prec", "fpr", "spec", "phi", "f")) {
  requireNamespace("dplyr")
  if (is.null(names(rankings))) {
    ranking.names <- seq_along(rankings)
  } else {
    ranking.names <- names(rankings)
  }
  evals <- lapply(rankings, evaluate.ranking, true.matrix=true.matrix, perf.measures=perf.measures)
  metrics <- dplyr::bind_rows(lapply(ranking.names, function(rn) {
    data.frame(ranking.name=rn, evals[[rn]]$metrics, check.names = F, stringsAsFactors = F)
  }))
  au.score <- dplyr::bind_rows(lapply(ranking.names, function(rn) {
    data.frame(ranking.name=rn, evals[[rn]]$au.score, check.names = F, stringsAsFactors = F)
  }))
  list(metrics=metrics, au.score=au.score)
}