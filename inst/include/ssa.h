#ifndef DYNGEN_SSA_H
#define DYNGEN_SSA_H

#include <Rcpp.h>

using namespace Rcpp;

class SSA {
public:
  SSA() {}
  
  void step(
      const NumericVector& state, 
      const NumericVector& transition_rates, 
      const NumericMatrix& nu,
      double* dtime, 
      NumericVector& dstate
  );
  
  void calculate_transition_rates(
      const NumericVector& state,
      const NumericVector& params,
      const double time,
      NumericVector& transition_rates
  );
  
  List simulate(
      const NumericVector& initial_state,
      const NumericVector& params,
      const NumericMatrix& nu,
      const double final_time,
      const double max_duration,
      const bool stop_on_neg_state,
      const bool stop_on_neg_propensity,
      const bool verbose,
      const double time_next_console
  );
} ;

#endif