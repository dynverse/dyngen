#' @importFrom stats setNames
generate_formulae <- function(net, geneinfo, cells=tibble(cell_id=1, dies=FALSE)) {
  formulae <- list()
  
  # generate the formulae for every target gene
  for (target_id in geneinfo$gene_id) {
    info <- dynutils::extract_row_to_list(geneinfo, which(geneinfo$gene_id == target_id))
    cell_id <- info$cell_id
    
    rg <- glue("rg_{cell_id}")
    kg <- glue("kg_{cell_id}")
    pg <- glue("pg_{cell_id}")
    qg <- glue("qg_{cell_id}")
    dg <- glue("dg_{cell_id}")
    
    y <- glue("y_{target_id}")
    q <- glue("q_{target_id}")
    
    # only add gene regulation if this is a gene
    if(info$isgene) {
      x <- glue("x_{target_id}")
      r <- glue("r_{target_id}")
      d <- glue("d_{target_id}")
      p <- glue("p_{target_id}")
      
      a0 <- glue("a0_{target_id}")
      
      regulator_ids <- net %>% filter(net$to == target_id) %>% pull(from)
      
      # regulator binding sites
      if (length(regulator_ids) > 0) {
        regulator_ys <- map(regulator_ids, ~glue("y_{.}"))
        regulator_ks <- map(regulator_ids, ~glue("k_{.}"))
        regulator_cs <- map(regulator_ids, ~glue("c_{.}"))
        
        # calculate the formulae for binding
        regulator_affinities <- pmap(list(regulator_ys, regulator_ks, regulator_cs), function(y, k, c) glue("(({y}/{k})**{c})"))
      
        # now combine the binding combinations
        binding_configurations <- map(get_configuration_ids(length(regulator_ids)), get_binding_configuration, length(regulator_ids))
        configuration_affinities <- map(binding_configurations, ~paste0(regulator_affinities[.], collapse="*"))
        configuration_as <- map(get_configuration_ids(length(regulator_ids)), ~glue("a_{target_id}_{.}"))
        
        activations <- map2(configuration_affinities, configuration_as, function(affinity, a) glue("{a} * {affinity}"))
        activation <- paste0(activations, collapse = "+")
        
        # now make the activation function, by summing the activations
        activation_function <- glue("({a0} + {activation})/(1 + {activation})")
      } else {
        activation_function <- "a0"
      }
      
      # mRNA production
      formula_id <- paste0("mrna_production_", target_id)
      formula <- glue("{rg} * {r} * {activation_function}")
      formulae[[formula_id]] <- list(
        formula_id = formula_id,
        formula = formula,
        effect = 1,
        molecule = x
      )
      
      # mRNA degradation
      formula_id <- paste0("mrna_degradation_", target_id)
      formula = glue("{dg} * {d} * {x}")
      formulae[[formula_id]] <- list(
        formula_id = formula_id,
        formula = formula,
        effect = -1,
        molecule = x
      )
    }
    
    # add protein synthesis if target gene also regulates other genes
    if(target_id %in% net$from) {
      if(info$isgene) {
        # protein production
        formula_id <- paste0("protein_production", target_id)
        formula = glue("{pg} * {p} * {x}")
        formulae[[formula_id]] <- list(
          formula_id = formula_id,
          formula = formula,
          effect = 1,
          molecule = y
        )
      }
      
      # protein degradation
      formula_id <- paste0("protein_degradation", target_id)
      formula = glue("{qg} * {q} * {y}")
      formulae[[formula_id]] <- list(
        formula_id = formula_id,
        formula = formula,
        effect = -1,
        molecule = y
      )
    }
  }
  
  # global transcription rates
  # for (i in 1:nrow(cells)) {
  #   cellinfo = cells[i,]
  #   
  #   rg = rgs[[celltypeinfo$celltype]]
  #   
  #   formula = fcon(paste0("ifelse(", rg@string, ">0.2, 1, 0) * 10"))
  #   nu = get_nu(higher=rg)
  #   add.formula(formula, nu, paste0("global_mrna_production_increase"), formulae)
  #   
  #   formula = fcon(paste0("1 * 10 * ", rg@string))
  #   nu = get_nu(lower=rg)
  #   add.formula(formula, nu, paste0("global_mrna_production_decrease"), formulae)
  #   
  #   if(celltypeinfo$dies) {
  #     formula = fcon(paste0(rg@string, " * (",x@string, "*100)^2"))
  #     nu = get_nu(lower=rg)
  #     add.formula(formula, nu, paste0("death"), formulae)
  #   }
  # }
  
  formulae %>% bind_rows()
}
