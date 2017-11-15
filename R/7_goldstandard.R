library(igraph)

smooth_expression <- function(expression, smooth_window=50) {
  expression %>% 
    zoo::rollmean(smooth_window, c("extend", "extend", "extend")) %>%
    magrittr::set_rownames(rownames(expression))
}

preprocess_simulation_for_gs <- function(simulation, model, smooth_window=50) {
  # print("smoothing...")
  
  geneinfo <- filter(model$geneinfo, main, !is.na(module_id))
  expression <- simulation$expression[, geneinfo %>% pull(gene_id)]
  
  simulation$expression_smooth <- expression %>% as.data.frame() %>% split(simulation$stepinfo$simulation_id) %>% map(smooth_expression, smooth_window=smooth_window) %>% do.call(rbind, .)
  dimnames(simulation$expression_smooth) <- dimnames(expression)
  
  # print("normalizing...")
  simulation$expression_normalized <- dynutils::scale_quantile(simulation$expression_smooth, outlier_cutoff=0.05)
  
  # print("calculating module expression...")
  simulation$expression_modules <- simulation$expression_normalized %>% t %>% as.data.frame() %>% split(factor(as.numeric(geneinfo$module_id), levels=model$modulenodes$module_id)) %>% map(~apply(., 2, mean)) %>% do.call(rbind, .) %>% t %>% magrittr::set_colnames(unique(geneinfo$module_id))
  simulation$expression_modules <- simulation$expression_modules[, as.character(model$modulenodes$module_id)] # fix ordering
  
  simulation
}

# load in milestone network and edge_operations
# model$edge_operations <- read_tsv(paste0("data/modulenetworks/", params$model$modulenet_name, "/edge_operations.tsv"), col_types=cols())

# extract all possible milestone paths from the start_milestone
# paths <- map(start_milestones, ~all_simple_paths(milestone_graph, as.character(.))) %>% unlist(recursive=FALSE) %>% map(~names(.) %>% as.numeric()) # this does not include cycles
get_milestone_paths <- function(milestone_network, start_milestones, max_path_length = 10) {
  recursor <- function(curnode, curpath=c(), max_path_length = 10) {
    curpath <- c(curpath, curnode)
    tos <- milestone_network %>% filter(from == curnode) %>% pull(to)
    
    if (length(curpath) < max_path_length) {
      c(list(curpath), unlist(map(tos, recursor, curpath=curpath), recursive=FALSE))
    } else {
      curpath
    }
  }
  paths <- map(start_milestones, recursor, max_path_length=max_path_length) %>% unlist(recursive=FALSE) %>% keep(~length(.) > 1)
  paths
}

process_operations <- function(edge_operations, module_ids) {
  edge_operations <- edge_operations %>% mutate(edge_id = seq_len(n()))
  
  # process operations
  operations <- edge_operations %>% 
    separate_rows(module_progression, sep="\\|") %>% 
    mutate(operation_id = seq_len(n())) %>% 
    separate_rows(module_progression, sep=",") %>% 
    mutate(
      operation = c(1, -1)[as.numeric(factor(substring(module_progression, 1, 1), levels=c("+", "-")))], 
      module_id = as.integer(substring(module_progression, 2))
    )
  operations$module_id <- factor(operations$module_id, levels=model$modulenodes$module_id) # factor -> acast
  operations
}

# now extract operations along every path
extract_path_operations <- function(operations, paths) {
  path_operations <- map(paths, function(path) {
    map2(path[-length(path)], path[-1], function(a, b) {
      operations %>% filter(from == a, to == b)
    }) %>% bind_rows()
  })
  path_operations
}

# extract reference expression for every edge
extract_references <- function(path_operations, milestone_network, reference_length = 100) {
  references <- map(path_operations, function(operations) {
    # put every consecutive operation id in a separate run
    operations$operation_run <- rle(operations$operation_id) %>% {rep(seq_along(.$length), .$length)}
    
    # first calculate the changes during every operation (run), than do a cumsum to get the actual reference
    reference_expression <- reshape2::acast(operations, module_id~operation_run, value.var="operation", fill=0, drop=FALSE) %>% apply(1, cumsum)
    
    if(class(reference_expression)== "numeric") {
      reference_expression <- matrix(reference_expression, nrow=1) %>% magrittr::set_colnames(names(reference_expression))
    }
    
    # add start
    reference_expression <- rbind(rep(0, ncol(reference_expression)), reference_expression)
    
    if (any(is.na(colnames(reference_expression)))) {stop("Check edge operations for wrong operations")}
    
    # check for one operation
    if(class(reference_expression)== "numeric") {
      stop("Need at least two operations")
    }
    
    # check for incorrect operations (leading to higher than 1 or lower than 0)
    if(min(reference_expression) < 0 || max(reference_expression) > 1) {
      print(operations$edge_id)
      print(reference_expression)
      stop("Inconsistent operations")
    }
    
    reference_expression <- reference_expression %>% apply(2, function(x) {
      reference_expression <- approx(1:length(x), x, seq(1, length(x), length.out=reference_length * length(x)))$y
    })
    
    operations <- bind_rows(operations %>% filter(row_number() == 1) %>% mutate(operation_id = 0), operations)
    reference_info <- operations %>% 
      select(from, to, edge_id) %>% 
      {.[rep(seq_len(nrow(.)), each=reference_length), ]} %>% 
      mutate(reference_row_id = seq_len(n())) %>% 
      group_by(edge_id) %>% 
      mutate(percentage = seq(0, 1, length.out=n())) %>% 
      ungroup()
    
    lst(
      reference_expression, 
      reference_info, 
      open_ended = (tail(operations$to, 1) %in% milestone_network$from) > 1
    )
  })
  references
}

# pheatmap::pheatmap(samplexpression, cluster_cols=F, cluster_rows=F)

#' @importFrom pdist pdist
#' @importFrom dtw dtw
map_to_reference <- function(simulation_expressions, references, ncores = 8) {
  times <- pbapply::pblapply(cl = ncores, simulation_expressions, function(simulation_expression) {
    dtws <- map(references, function(reference) {
      if(ncol(simulation_expression) != ncol(reference$reference_expression)) {
        print(colnames(simulation_expression))
        print(colnames(reference$reference_expression))
        stop("Colnames of reference and simulation should be equal!")
      }
      
      simulation_expression <- simulation_expression[, colnames(reference$reference_expression)]
      
      dist <- pdist::pdist(simulation_expression, reference$reference_expression) %>% as.matrix()
      dtw <- dtw::dtw(dist, open.end = TRUE)
      dtw
    })
    
    normalizedDistance <- min(which.min(map_dbl(dtws, "normalizedDistance")))
    
    best_reference_id <- which.min(map_dbl(dtws, "normalizedDistance"))
    dtw <- dtws[[best_reference_id]]
    reference <- references[[best_reference_id]]
    
    # if open ended, add extra unassigned index1 to last index2 
    dtw$index1 <- c(dtw$index1, which(!seq_len(nrow(simulation_expression)) %in% dtw$index1))
    dtw$index2 <- c(dtw$index2, rep(tail(dtw$index2, 1), length(dtw$index2) - length(dtw$index1)))
    
    times <- tibble(index1=dtw$index1, index2=dtw$index2) %>% 
      mutate(step_id = rownames(simulation_expression)[index1]) %>% 
      group_by(step_id) %>% 
      summarise(reference_row_id = round(mean(index2))) %>% 
      left_join(reference$reference_info, by="reference_row_id") %>% 
      select(-reference_row_id) %>% 
      mutate(percentage = ifelse(is.na(percentage), 1, percentage)) %>% 
      mutate(normalizedDistance = normalizedDistance)
    times
  }) %>% bind_rows()
  times
}