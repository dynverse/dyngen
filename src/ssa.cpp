#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

List resize( const List& x, int n ){
  int oldsize = x.size() ;
  List y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}

List SSA::simulate(
    const NumericVector& initial_state,
    const NumericVector& params,
    const NumericMatrix& nu,
    const double final_time,
    const double max_duration,
    const bool stop_on_neg_state,
    const bool stop_on_neg_propensity,
    const bool verbose,
    const double time_next_console
) {
  List output(10);
  
  NumericVector state = clone(initial_state);
  double time = 0;
  NumericVector transition_rates = no_init(nu.ncol());
  
  // calculate initial transition rates
  calculate_transition_rates(state, params, time, transition_rates);
  
  output(0) = List::create(
    _["time"] = time,
    _["state"] = clone(state),
    _["transition_rates"] = clone(transition_rates)
  );
  
  return output;
}

void SSA::calculate_transition_rates(
    const NumericVector& state,
    const NumericVector& params,
    const double time,
    NumericVector& transition_rates
) {
  // this function should be overridden
}
