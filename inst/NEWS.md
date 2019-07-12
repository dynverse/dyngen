# dyngen 0.2.0 (2019-07-12)

 * Complete rewrite of `dyngen`. Major improvements: 
 
   - All aspects of the pipeline have been optimised towards execution time and end-user usability.
   - `dyngen` 0.2.0 uses `fastgssa` 0.2.0, which has also been rewritten entirely in `Rcpp`,
     thereby improving the speed significantly.
   - Mapping of a simulation to the gold standard is greatly improved.
   - Custom backbones can be defined using backbone lego pieces. See `?bblego` for more information.

# dyngen 0.1.0 (2017-04-27)

 * Initial release of `dyngen`, a package for generating synthetic single-cell data from regulatory networks.
   Key features are:
   
   - The cells undergo a dynamic process throughout the simulation.
   - Many different trajectory types are supported.
   - `dyngen` 0.1.0 uses `fastgssa` 0.1.0, a clone of `GillespieSSA` that is much less
     error-prone and more efficient than `GillespieSSA`.

# dyngen 0.0.1 (2016-04-04)

 * Just a bunch of scripts on a repository, which creates random networks using `igraph` and 
   generates simple single-cell expression data using `GillespieSSA`.