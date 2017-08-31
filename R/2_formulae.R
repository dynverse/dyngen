add.formula = function(formula, nu, name, system) {system$formulae[[name]] <- formula;system$nus.changes[[name]] <- nu;assign("system", system, parent.frame())}

add.variable = function(variable, ..., system) {
  varname = variable@string
  info = list(name=varname, ...)
  system$variables[[varname]] <- info
  
  group = str_replace(varname, "^(.*?)[_$].*", "\\1")
  if (group %in% c("b", "u", "x", "y", "k", "rg", "r", "d", "p", "q", "a0", "a", "assoc", "deassoc")) {
    system$vargroups[[group]] <- c(system$vargroups[[group]], varname)
  }
  
  assign("system", system, parent.frame())
  variable
}

generate_formulae = function(net, genes, celltypes=tibble(celltype=1, dies=F)) {
  ## generating reaction formulas (the probability that a reaction occurs in [t, t+dt])
  
  system = list(formulae=list(), nus.changes=list(), variables=list(), vargroups=list())
  
  get.nu = function(higher=NULL, lower=NULL) {setNames(c(rep(1, length(higher)), rep(-1, length(lower))), lapply(c(higher, lower), function(x) {x@string}))}
  
  kg = add.variable(fvar("kg"), system=system) # global tf binding speed (influences burstiness of transcriptional regulation)
  rgs = c()
  for(celltype in celltypes$celltype) {
    rgs = c(rgs, add.variable(fvar(paste0("rg_", celltype)), system=system))# global transcription rate
  }
  dg = add.variable(fvar("dg"), system=system) # global RNA degradation rate
  pg = add.variable(fvar("pg"), system=system) # global translation rate
  qg = add.variable(fvar("qg"), system=system) # global protein degradation rate
  
  for (target in genes$gene) {
    geneinfo = genes[genes$gene == target,] %>% as.list
    
    y = add.variable(fvar("y", target), gene=target, type="protein_production", system=system)
    q = add.variable(fvar("q", target), gene=target, type="protein_degradation", system=system)
    
    # only add gene regulation if this is a gene
    if(geneinfo$isgene) {
      x = add.variable(fvar("x", target), gene=target, system=system)
      r = add.variable(fvar("r", target), gene=target, system=system)
      d = add.variable(fvar("d", target), gene=target, system=system)
      
      p = add.variable(fvar("p", target), gene=target, system=system)
      a0 = add.variable(fvar("a0", target), gene=target, system=system)
      
      subnet = net[net$to == target,]
      regs = subnet$from
      effects = subnet$effect
      
      if(!all(effects %in% c(1, -1))) stop("effects should be either 1 or -1")
      
      rg = rgs[[geneinfo$celltype]]
      
      boundvars = list()
      
      if(!is.na(geneinfo$a0)) {
        a0_value = geneinfo$a0
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
          regulator = regs[[i]]
          effect = effects[[i]]
          strength = subnet[i,]$strength
          cooperativity = ifelse(!is.na(subnet[i,]$cooperativity), subnet[i,]$cooperativity, 2)
          
          #u = add.variable(fvar("u", target, regulator), target=target, regulator=regulator)
          #b = add.variable(fvar("b", target, regulator), target=target, regulator=regulator)
          k = add.variable(fvar("k", target, regulator), target=target, regulator=regulator, strength=strength, system=system, randomize=subnet[i,]$randomize)
          
          #boundvars = c(boundvars, b)
          
          y_regulator = fvar("y", regulator)
          
          # binding
          #if (cooperativity > 1) {
            #formula = fcon(paste0(u@string, " * ", y_regulator@string, "^", cooperativity, " * ", kg@string))
            #formula = u * y_regulator * y_regulator * kg
          #} else {
            #formula = u * y_regulator * kg
          #}
          
          #nu = get.nu(b, list(y_regulator, u))
          #add.formula(formula, nu)
          
          # deassociation
          #formula = b * k * kg
          #nu = get.nu(list(y_regulator, u), b)
          #add.formula(formula, nu)
          
          # add to decision tree
          a = add.variable(fvar("a", target, regulator), target=target, regulator=regulator, effect=effect, system=system)
          
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
      formula = rg * r * decisiontree
      nu = get.nu(x)
      add.formula(formula, nu, paste0("mrna_production_", target), system)
      
      # mRNA degradation
      formula = dg * d * x
      nu = get.nu(lower=x)
      add.formula(formula, nu, paste0("mrna_degradation_", target), system)
    }
    
    if(target %in% net$from || TRUE) {
      if(geneinfo$isgene) {
        # protein production
        formula = pg * p * x
        nu = get.nu(y)
        add.formula(formula, nu, paste0("protein_production_", target), system)
      }
      
      # protein degradation
      formula = qg * q * y
      nu = get.nu(lower=y)
      add.formula(formula, nu, paste0("protein_degradation_", target), system)
    }
  }
  
  for (i in 1:nrow(celltypes)) {
    celltypeinfo = celltypes[i,]
    
    rg = rgs[[celltypeinfo$celltype]]
    
    formula = fcon(paste0("ifelse(", rg@string, ">0.2, 1, 0) * 10"))
    nu = get.nu(higher=rg)
    add.formula(formula, nu, paste0("global_mrna_production_increase"), system)
    
    formula = fcon(paste0("1 * 10 * ", rg@string))
    nu = get.nu(lower=rg)
    add.formula(formula, nu, paste0("global_mrna_production_decrease"), system)
    
    if(celltypeinfo$dies) {
      formula = fcon(paste0(rg@string, " * (",x@string, "*100)^2"))
      nu = get.nu(lower=rg)
      add.formula(formula, nu, paste0("death"), system)
    }
  }
  
  system$formulae.strings <- map_chr(system$formulae, function(fl) fl@string)
  
  system
}
