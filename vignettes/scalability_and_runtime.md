On runtime and scalability
================

``` r
library(tidyverse)
library(dyngen)
```

<!-- github markdown built using 
rmarkdown::render("vignettes/advanced_constructing_backbone.Rmd", output_format = rmarkdown::github_document())
-->

In this vignette, we will take a look at the runtime of dyngen as the
number of genes and the number of cells sampled is increased. To this
end, we’ll be using the bifurcating cycle backbone which is well known
for its beautiful 3D butterfly shape!

An example of a resulting dyngen model is shown here. We tweaked some of
the parameters by running this particular backbone once with
`num_cells = 100` and `num_features = 100` and verifying that the new
parameters still yield the desired outcome. The parameters we tweaked
are:

-   On average, 10 cells are sampled per simulation
    (`num_simulations = 100` and `num_cells = 1000`). You could increase
    this ratio to get a better cell count yield from a given set of
    simulations, but cells from the same simulation that are temporally
    close will have highly correlated expression profiles.
-   Increased time steps `tau`. This will make the Gillespie algorithm
    slighty faster but might result in unexpected artifacts in the
    simulated data.
-   `census_interval` increased from 4 to 10. This will cause dyngen to
    store an expression profile only every 10 time units. Since the
    total simulation time is xxx, each simulation will result in yyy
    data points. Note that on average only 10 data points are sampled
    per simulation.

## Experiments

``` r
library(dyngen)
library(tidyverse)

set.seed(1)

save_dir <- "scalability_and_runtime_runs"
dir.create(save_dir, recursive = TRUE)
```

    ## Warning in dir.create(save_dir, recursive = TRUE): 'scalability_and_runtime_runs' already exists

``` r
backbone <- backbone_bifurcating_cycle()

settings <- bind_rows(
  tibble(num_cells = 10000, num_features = 10000, rep = 1),
  crossing(
    num_cells = seq(1000, 10000, by = 1000),
    num_features = 100,
    rep = seq_len(3)
  ),
  crossing(
    num_cells = 100,
    num_features = seq(1000, 10000, by = 1000),
    rep = seq_len(3)
  )
) %>% 
  mutate(filename = paste0(save_dir, "/cells", num_cells, "_feats", num_features, "_rep", rep, ".rds"))

timings <- pmap_dfr(settings, function(num_cells, num_features, rep, filename) {
  if (!file.exists(filename)) {
    cat("Running num_cells: ", num_cells, ", num_features: ", num_features, ", rep: ", rep, "\n", sep = "")
    num_tfs <- nrow(backbone$module_info)
    num_targets <- round((num_features - num_tfs) / 2)
    num_hks <- num_features - num_targets - num_tfs
    
    out <- 
      initialise_model(
        backbone = backbone,
        num_tfs = num_tfs,
        num_targets = num_targets,
        num_hks = num_hks,
        num_cells = num_cells,
        gold_standard_params = gold_standard_default(
          census_interval = 1,
          tau = 100/3600
        ),
        simulation_params = simulation_default(
          census_interval = 10,
          ssa_algorithm = ssa_etl(tau = 300/3600),
          experiment_params = simulation_type_wild_type(
            num_simulations = num_cells / 10
          )
        ),
        verbose = FALSE
      ) %>% 
      generate_dataset()
    
    tim <- 
      get_timings(out$model) %>% 
      mutate(rep, num_cells, num_features)
    
    write_rds(tim, filename, compress = "gz")
  }
  
  read_rds(filename)
})
```

## Single run, 10k cells, 10k features

``` r
timings0 <- 
  timings %>% 
  filter(num_cells == 10000, num_features == 10000) %>% 
  mutate(name = forcats::fct_rev(forcats::fct_inorder(paste0(group, ": ", task))))

ggplot(timings0) + 
  geom_bar(aes(x = name, y = time_elapsed, fill = group), stat = "identity") +
  scale_fill_brewer(palette = "Dark2") + 
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip() + 
  labs(x = NULL, y = "Time (s)", fill = "dyngen stage")
```

![](scalability_and_runtime_files/figure-gfm/bblego-1.png)<!-- -->

## Increasing the number of cells

``` r
timings1 <- 
  timings %>% 
  filter(num_features == 100)

timings1 %>% 
  group_by(group, task, num_cells, num_features) %>% 
  summarise(time_elapsed = mean(time_elapsed), groups = "drop") %>% 
  group_by(num_cells, num_features, group) %>% summarise(time_elapsed = sum(time_elapsed), .groups = "drop") %>% 
  ggplot() + 
  geom_bar(aes(x = forcats::fct_inorder(as.character(num_cells)), y = time_elapsed, fill = forcats::fct_inorder(group)), stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Number of cells", y = "Average time (s)", fill = "dyngen step")
```

    ## `summarise()` regrouping output by 'group', 'task', 'num_cells' (override with `.groups` argument)

![](scalability_and_runtime_files/figure-gfm/figure1-1.png)<!-- -->

## Increasing the number of features

``` r
timings2 <- 
  timings %>% 
  filter(num_cells == 100)

timings2 %>% 
  group_by(group, task, num_cells, num_features) %>% 
  summarise(time_elapsed = mean(time_elapsed), groups = "drop") %>% 
  group_by(num_cells, num_features, group) %>% summarise(time_elapsed = sum(time_elapsed), .groups = "drop") %>% 
  ggplot() + 
  geom_bar(aes(x = forcats::fct_inorder(as.character(num_features)), y = time_elapsed, fill = forcats::fct_inorder(group)), stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Number of features", y = "Average time (s)", fill = "dyngen step")
```

    ## `summarise()` regrouping output by 'group', 'task', 'num_cells' (override with `.groups` argument)

![](scalability_and_runtime_files/figure-gfm/figure2-1.png)<!-- -->

## Ratio between `num_cells` and `num_simulations` is too large

## `tau` is too large

## Execution platform

These timings were measured using 30 (out of 32) threads using a AMD
Ryzen 9 5950X clocked at 3.4GHz.

Session info:

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-redhat-linux-gnu (64-bit)
    ## Running under: Fedora 33 (Workstation Edition)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib64/libflexiblas.so.3.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] forcats_0.5.0   stringr_1.4.0   dplyr_1.0.2     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2     tibble_3.0.5    ggplot2_3.3.3   tidyverse_1.3.0 dyngen_0.4.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fs_1.5.0            dynutils_1.0.5.9000 lubridate_1.7.9.2   RColorBrewer_1.1-2  httr_1.4.2          tools_4.0.3         backports_1.2.1     R6_2.5.0            irlba_2.3.3         DBI_1.1.1           colorspace_2.0-0    withr_2.4.0         tidyselect_1.1.0    carrier_0.1.0      
    ## [15] gridExtra_2.3       processx_3.4.5      proxyC_0.1.5        compiler_4.0.3      cli_2.2.0           rvest_0.3.6         RcppXPtrUtils_0.1.1 xml2_1.3.2          labeling_0.4.2      scales_1.1.1        pbapply_1.4-3       digest_0.6.27       rmarkdown_2.6       pkgconfig_2.0.3    
    ## [29] htmltools_0.5.1     dbplyr_2.0.0        rlang_0.4.10        readxl_1.3.1        rstudioapi_0.13     farver_2.0.3        generics_0.1.0      jsonlite_1.7.2      magrittr_2.0.1      GillespieSSA2_0.2.7 patchwork_1.1.1     Matrix_1.2-18       Rcpp_1.0.6          munsell_0.5.0      
    ## [43] fansi_0.4.2         viridis_0.5.1       lifecycle_0.2.0     stringi_1.5.3       yaml_2.2.1          ggraph_2.0.4        MASS_7.3-53         plyr_1.8.6          grid_4.0.3          parallel_4.0.3      ggrepel_0.9.0       crayon_1.3.4.9000   lattice_0.20-41     dynwrap_1.2.2      
    ## [57] graphlayouts_0.7.1  haven_2.3.1         hms_1.0.0           knitr_1.30          ps_1.5.0            pillar_1.4.7        igraph_1.2.6        babelwhale_1.0.1    reshape2_1.4.4      codetools_0.2-16    reprex_0.3.0        glue_1.4.2          evaluate_0.14       remotes_2.2.0      
    ## [71] RcppParallel_5.0.2  modelr_0.1.8        vctrs_0.3.6         dynparam_1.0.1      tweenr_1.0.1        cellranger_1.1.0    gtable_0.3.0        polyclip_1.10-0     assertthat_0.2.1    xfun_0.20           ggforce_0.3.2       lmds_0.1.0          broom_0.7.2         tidygraph_1.2.0    
    ## [85] viridisLite_0.3.0   ellipsis_0.3.1
