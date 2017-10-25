add.formula = function(formula, nu, name, system) {system$formulae[[name]] <- formula;system$nus <- c(system$nus, list(nu));assign("system", system, parent.frame())}

add.variable = function(variable, ..., system) {
  varname <- variable@string
  info <- list(name=varname, ...)
  
  group <- str_replace(varname, "^(.*?)[_$].*", "\\1")
  system$vargroups[[group]] <- c(system$vargroups[[group]], varname)
  info$group <- group
  
  system$variables[[varname]] <- info
  
  assign("system", system, parent.frame())
  variable
}

#' @importFrom stats setNames
generate_formulae <- function(net, geneinfo, cells=tibble(cell_id=1, dies=FALSE)) {
  # initialize system
  system <- list(formulae=list(), nus=list(), variables=list(), vargroups=list())
  
  # helper function to quickly get the nu based on synthesis and degradation
  get_nu <- function(higher=NULL, lower=NULL) {
    bind_rows(
      tibble(effect = 1, formula_id=map_chr(higher, ~.@string)),
      tibble(effect = -1, formula_id=map_chr(lower, ~.@string))
    )
  }
  
  # variables for every cell in the system
  kgs <- c()
  rgs <- c()
  pgs <- c()
  qgs <- c()
  dgs <- c()
  for(cell_id in cells$cell_id) {
    kgs <- c(kgs, add.variable(fvar(paste0("kg_", cell_id)), system=system, cell_id=cell_id)) # global tf binding
    rgs = c(rgs, add.variable(fvar(paste0("rg_", cell_id)), system=system, cell_id=cell_id)) # global transcription rate
    dgs = c(dgs, add.variable(fvar(paste0("dg_", cell_id)), system=system, cell_id=cell_id)) # global mRNA degradation rate
    pgs = c(pgs, add.variable(fvar(paste0("pg_", cell_id)), system=system, cell_id=cell_id)) # global protein production rate
    qgs = c(qgs, add.variable(fvar(paste0("qg_", cell_id)), system=system, cell_id=cell_id)) # global protein degradation rate
  }
  
  # generate the formulae for every target gene
  for (target_id in geneinfo$gene_id) {
    info <- geneinfo[geneinfo$gene_id == target_id, ] %>% as.list()
    
    cell_id <- info$cell_id
    rg = rgs[[info$cell_id]]
    kg = kgs[[info$cell_id]]
    pg = pgs[[info$cell_id]]
    qg = qgs[[info$cell_id]]
    dg = dgs[[info$cell_id]]
    
    y = add.variable(fvar("y", target_id), gene_id=target_id, system=system)
    q = add.variable(fvar("q", target_id), gene_id=target_id, system=system)
    
    # only add gene regulation if this is a gene
    if(info$isgene) {
      x = add.variable(fvar("x", target_id), gene=target_id, system=system)
      r = add.variable(fvar("r", target_id), gene=target_id, system=system)
      d = add.variable(fvar("d", target_id), gene=target_id, system=system)
      
      p = add.variable(fvar("p", target_id), gene=target_id, system=system)
      a0 = add.variable(fvar("a0", target_id), gene=target_id, system=system)
      
      subnet = net[net$to == target_id,]
      regs = subnet$from
      effects = subnet$effect
      
      boundvars = list()
      
      if(!is.na(info$a0)) {
        a0_value = info$a0
      } else if(length(regs) == 0) {
        a0_value = 1
      } else if(sum(effects == -1) == 0) {
        a0_value = 0.0001
      } else if(sum(effects == 1) == 0) {
        a0_value = 1
      } else {
        a0_value = 0.5
      }
      system$variables[[a0@string]]$value = a0_value
      
      # regulator binding sites
      if (length(regs) > 0) {
        for (i in seq_len(length(regs))) {
          regulator_id = regs[[i]]
          effect = effects[[i]]
          strength = subnet[i,]$strength
          cooperativity = ifelse(!is.na(subnet[i,]$cooperativity), subnet[i,]$cooperativity, 2)
          
          #u = add.variable(fvar("u", target, regulator), target=target, regulator=regulator)
          #b = add.variable(fvar("b", target, regulator), target=target, regulator=regulator)
          k = add.variable(fvar("k", target_id, regulator_id), target_id=target_id, regulator_id=regulator_id, strength=strength, system=system, randomize=subnet[i,]$randomize)
          
          #boundvars = c(boundvars, b)
          
          y_regulator = fvar("y", regulator_id)
          
          # binding
          #if (cooperativity > 1) {
            #formula = fcon(paste0(u@string, " * ", y_regulator@string, "^", cooperativity, " * ", kg@string))
            #formula = u * y_regulator * y_regulator * kg
          #} else {
            #formula = u * y_regulator * kg
          #}
          
          #nu = get_nu(b, list(y_regulator, u))
          #add.formula(formula, nu)
          
          # deassociation
          #formula = b * k * kg
          #nu = get_nu(list(y_regulator, u), b)
          #add.formula(formula, nu)
          
          # add to decision tree
          a = add.variable(fvar("a", target_id, regulator_id), target=target_id, regulator=regulator_id, effect=effect, system=system)
          
          # quick and dirty decision tree, the last interaction will have priority
          if(i == 1) {
            decisiontree = fcon(paste0(a0@string, "+", a@string, " * (", y_regulator@string, "/", k@string, ")^", cooperativity))
            decisiontree_denom = fcon(paste0(1, "+ (", y_regulator@string, "/", k@string, ")^", cooperativity))
            
            #decisiontree = fcon(paste0(a0@string, "+", b@string, "/", u@string, "*", a@string))
            #decisiontree_denom = fcon(paste0(1, "+", b@string, "/", u@string))
            
            #decisiontree = fcon(paste0(a0@string, "+", b@string, "*", a@string))
            #decisiontree_denom = fcon(paste0(1, "+", b@string))
          } else {
            decisiontree = fcon(paste0(decisiontree@string, "+", a@string, " * (", y_regulator@string, "/", k@string, ")^", cooperativity))
            decisiontree_denom = fcon(paste0(decisiontree_denom@string, "+ (", y_regulator@string, "/", k@string, ")^", cooperativity))
            
            #decisiontree = fcon(paste0(decisiontree@string, "+", b@string, "/", u@string, "*", a@string))
            #decisiontree_denom = fcon(paste0(decisiontree_denom@string, "+", b@string, "/", u@string))
            
            #decisiontree = fcon(paste0(decisiontree@string, "+", b@string, "*", a@string))
            #decisiontree_denom = fcon(paste0(decisiontree_denom@string, "+", b@string))
          }
        }
      } else {
        decisiontree = a0
        decisiontree_denom = fcon(1)
      }
      
      decisiontree = fcon(paste0("(", decisiontree@string, ")/(", decisiontree_denom@string, ")"))
      
      # mRNA production
      formula_id <- paste0("mrna_production_", target_id)
      formula = rg * r * decisiontree
      nu = tibble(formula_id=formula_id, effect=1, molecule=x@string)
      add.formula(formula, nu, formula_id, system)
      
      # mRNA degradation
      formula_id <- paste0("mrna_degradation_", target_id)
      formula = dg * d * x
      nu = tibble(formula_id=formula_id, effect=-1, molecule=x@string)
      add.formula(formula, nu, formula_id, system)
    }
    
    # add protein synthesis if target gene also regulates other genes
    if(target_id %in% net$from) {
      if(info$isgene) {
        # protein production
        formula_id <- paste0("protein_production_", target_id)
        formula = pg * p * x
        nu = tibble(formula_id=formula_id, effect=1, molecule=y@string)
        add.formula(formula, nu, formula_id, system)
      }
      
      # protein degradation
      formula_id <- paste0("protein_degradation_", target_id)
      formula = qg * p * y
      nu = tibble(formula_id=formula_id, effect=-1, molecule=y@string)
      add.formula(formula, nu, formula_id, system)
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
  #   add.formula(formula, nu, paste0("global_mrna_production_increase"), system)
  #   
  #   formula = fcon(paste0("1 * 10 * ", rg@string))
  #   nu = get_nu(lower=rg)
  #   add.formula(formula, nu, paste0("global_mrna_production_decrease"), system)
  #   
  #   if(celltypeinfo$dies) {
  #     formula = fcon(paste0(rg@string, " * (",x@string, "*100)^2"))
  #     nu = get_nu(lower=rg)
  #     add.formula(formula, nu, paste0("death"), system)
  #   }
  # }
  
  system$formulae_strings <- map_chr(system$formulae, function(fl) fl@string)
  system$nus_df <- bind_rows(system$nus)
  
  system[names(system) != "nus"]
}
