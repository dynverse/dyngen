Advanced: On scalability and runtime
================

``` r
library(tidyverse)
library(dyngen)
```

<!-- github markdown built using 
rmarkdown::render("vignettes/advanced_scalability_and_runtime.Rmd", output_format = rmarkdown::github_document())
-->

In this vignette, we will take a look at the runtime of dyngen as the
number of genes and the number of cells sampled is increased. We’ll be
using the bifurcating cycle backbone which is well known for its
beautiful 3D butterfly shape!

``` r
library(dyngen)
library(tidyverse)

set.seed(1)

save_dir <- "advanced_scalability_and_runtime_runs"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

backbone <- backbone_bifurcating_cycle()
```

## Initial run

We’ll be running this simulation a few times, with different values for
`num_cells` and `num_features` to assess the scalability of dyngen. An
example of a resulting dyngen model is shown here.

``` r
num_cells <- 100
num_features <- 100
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
  generate_dataset(make_plots = TRUE)
```

    ## Loading required namespace: dynwrap

``` r
out$plot
```

![](advanced_scalability_and_runtime_files/figure-gfm/example-1.png)<!-- -->

We tweaked some of the parameters by running this particular backbone
once with `num_cells = 100` and `num_features = 100` and verifying that
the new parameters still yield the desired outcome. The parameters we
tweaked are:

-   On average, 10 cells are sampled per simulation
    (e.g. `num_simulations = 100` and `num_cells = 1000`). You could
    increase this ratio to get a better cell count yield from a given
    set of simulations, but cells from the same simulation that are
    temporally close will have highly correlated expression profiles.
-   Increased time steps `tau`. This will make the Gillespie algorithm
    slightly faster but might result in unexpected artifacts in the
    simulated data.
-   `census_interval` increased from 4 to 10. This will cause dyngen to
    store an expression profile only every 10 time units. Since the
    total simulation time is xxx, each simulation will result in yyy
    data points. Note that on average only 10 data points are sampled
    per simulation.

For more information on parameter tuning, see the vignette ‘Advanced:
tuning the simulation parameters’.

## Timing experiments

The simulations are run once with a large `num_features` and
`num_cells`, a few times with varying `num_cells` and then once more
with varying `num_features`. Every run is repeated three times in order
to get a bit more stable time measurements. Since some of the
simulations can take over 10 minutes, the timings results of the
simulations are cached in the
‘advanced\_scalability\_and\_runtime\_runs’ folder.\`

``` r
settings <- bind_rows(
  tibble(num_cells = 10000, num_features = 10000, rep = 1), #, rep = seq_len(3)),
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
    set.seed(rep)
    
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

timings_gr <- 
  timings %>% 
  group_by(group, task, num_cells, num_features) %>% 
  summarise(time_elapsed = mean(time_elapsed), .groups = "drop")

timings_sum <-
  timings %>% 
  group_by(num_cells, num_features, rep) %>%
  summarise(time_elapsed = sum(time_elapsed), .groups = "drop")
```

## Simulate a large dataset (10k × 10k)

Below is shown the timings of each of the steps in simulating a dyngen
dataset containing 10’000 genes and 10’000 features. The total
simulation time required is 1147 seconds, most of which is spent
performing the simulations itself.

``` r
timings0 <- 
  timings_gr %>% 
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

![](advanced_scalability_and_runtime_files/figure-gfm/bblego-1.png)<!-- -->

## Increasing the number of cells

By increasing the number of cells from 1000 to 10’000 whilst keeping the
number of features fixed, we can get an idea of how the simulation time
scales w.r.t. the number of cells.

``` r
timings1 <- 
  timings_gr %>% 
  filter(num_features == 100) %>% 
  group_by(num_cells, num_features, group) %>%
  summarise(time_elapsed = sum(time_elapsed), .groups = "drop")

ggplot(timings1) + 
  geom_bar(aes(x = forcats::fct_inorder(as.character(num_cells)), y = time_elapsed, fill = forcats::fct_inorder(group)), stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Number of cells", y = "Average time (s)", fill = "dyngen step")
```

![](advanced_scalability_and_runtime_files/figure-gfm/figure1-1.png)<!-- -->

It seems the execution time scales linearly w.r.t. the number of cells.
This makes sense, because as the number of cells are increased, so do we
increase the number of simulations made (which is not necessarily
mandatory). Since the simulations are independent of each other and take
up the most time, the execution time will scale linearly.

``` r
ggplot(timings_sum %>% filter(num_features == 100)) + 
  theme_bw() +
  geom_point(aes(num_cells, time_elapsed)) +
  scale_x_continuous(limits = c(0, 10000)) +
  scale_y_continuous(limits = c(0, 300)) +
  geom_abline(intercept = 22.097, slope = 0.0252) +
  labs(x = "Number of cells", y = "Execution time (s)")
```

![](advanced_scalability_and_runtime_files/figure-gfm/plot_timings_cell-1.png)<!-- -->

## Increasing the number of features

By increasing the number of features from 1000 to 10’000 whilst keeping
the number of cells fixed, we can get an idea of how the simulation time
scales w.r.t. the number of features

``` r
timings2 <- 
  timings_gr %>% 
  filter(num_cells == 100) %>% 
  group_by(num_cells, num_features, group) %>% 
  summarise(time_elapsed = sum(time_elapsed), .groups = "drop")

ggplot(timings2) + 
  geom_bar(aes(x = forcats::fct_inorder(as.character(num_features)), y = time_elapsed, fill = forcats::fct_inorder(group)), stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Number of features", y = "Average time (s)", fill = "dyngen step")
```

![](advanced_scalability_and_runtime_files/figure-gfm/figure2-1.png)<!-- -->

It seems the execution time also scales linearly w.r.t. the number of
features. As more genes are added to the underlying gene regulatory
network, the density of the graph doesn’t change, so it makes sense that
the execution time also scales linearly w.r.t. the number of features.

``` r
ggplot(timings_sum %>% filter(num_cells == 100)) + 
  theme_bw() +
  geom_point(aes(num_features, time_elapsed)) +
  scale_x_continuous(limits = c(0, 10000)) +
  scale_y_continuous(limits = c(0, 850)) +
  geom_abline(intercept = 0.5481, slope = 0.07988) +
  labs(x = "Number of features", y = "Execution time (s)")
```

![](advanced_scalability_and_runtime_files/figure-gfm/plot_timings_feats-1.png)<!-- -->

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
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8      
    ##  [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] dyngen_0.4.1    forcats_0.5.0   stringr_1.4.0   dplyr_1.0.2     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2     tibble_3.0.5    ggplot2_3.3.3   tidyverse_1.3.0
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1                backports_1.2.1             babelwhale_1.0.1            plyr_1.8.6                  igraph_1.2.6                proxyC_0.1.5                splines_4.0.3              
    ##   [8] dynwrap_1.2.2               BiocParallel_1.24.1         GenomeInfoDb_1.26.2         digest_0.6.27               htmltools_0.5.1             viridis_0.5.1               fansi_0.4.2                
    ##  [15] magrittr_2.0.1              memoise_1.1.0               carrier_0.1.0               limma_3.44.3                remotes_2.2.0               annotate_1.68.0             graphlayouts_0.7.1         
    ##  [22] modelr_0.1.8                RcppParallel_5.0.2          matrixStats_0.57.0          dynutils_1.0.5.9000         colorspace_2.0-0            blob_1.2.1                  rvest_0.3.6                
    ##  [29] ggrepel_0.9.0               haven_2.3.1                 xfun_0.20                   crayon_1.3.4.9000           RCurl_1.98-1.2              jsonlite_1.7.2              genefilter_1.72.1          
    ##  [36] survival_3.2-7              glue_1.4.2                  countsimQC_1.8.0            polyclip_1.10-0             gtable_0.3.0                GillespieSSA2_0.2.7         zlibbioc_1.36.0            
    ##  [43] XVector_0.30.0              DelayedArray_0.16.0         BiocGenerics_0.36.0         dynparam_1.0.1              scales_1.1.1                DBI_1.1.1                   edgeR_3.30.3               
    ##  [50] Rcpp_1.0.6                  viridisLite_0.3.0           xtable_1.8-4                bit_4.0.4                   RcppXPtrUtils_0.1.1         stats4_4.0.3                randtests_1.0              
    ##  [57] DT_0.16                     htmlwidgets_1.5.2           httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.1              pkgconfig_2.0.3             XML_3.99-0.5               
    ##  [64] farver_2.0.3                dbplyr_2.0.0                locfit_1.5-9.4              labeling_0.4.2              tidyselect_1.1.0            rlang_0.4.10                reshape2_1.4.4             
    ##  [71] AnnotationDbi_1.50.3        munsell_0.5.0               cellranger_1.1.0            tools_4.0.3                 cli_2.2.0                   generics_0.1.0              RSQLite_2.2.2              
    ##  [78] broom_0.7.2                 evaluate_0.14               yaml_2.2.1                  processx_3.4.5              knitr_1.30                  bit64_4.0.5                 fs_1.5.0                   
    ##  [85] tidygraph_1.2.0             lmds_0.1.0                  caTools_1.18.1              ggraph_2.0.4                pbapply_1.4-3               xml2_1.3.2                  compiler_4.0.3             
    ##  [92] rstudioapi_0.13             reprex_0.3.0                tweenr_1.0.1                geneplotter_1.68.0          stringi_1.5.3               ps_1.5.0                    lattice_0.20-41            
    ##  [99] Matrix_1.2-18               vctrs_0.3.6                 pillar_1.4.7                lifecycle_0.2.0             bitops_1.0-6                irlba_2.3.3                 patchwork_1.1.1            
    ## [106] GenomicRanges_1.42.0        R6_2.5.0                    gridExtra_2.3               IRanges_2.24.0              codetools_0.2-16            MASS_7.3-53                 assertthat_0.2.1           
    ## [113] SummarizedExperiment_1.20.0 DESeq2_1.30.0               rprojroot_2.0.2             withr_2.4.0                 S4Vectors_0.28.0            GenomeInfoDbData_1.2.4      parallel_4.0.3             
    ## [120] hms_1.0.0                   grid_4.0.3                  rmarkdown_2.6               MatrixGenerics_1.2.0        ggforce_0.3.2               Biobase_2.50.0              lubridate_1.7.9.2
