#include <Rcpp.h>
#include "ssa_direct.h"

using namespace Rcpp;

int weighted_sample(const NumericVector& weight) {
  double max = sum(weight);
  double ran = runif(1, 0, max)(0);
  int j = 0;
  while (ran > weight(j)) {
    ran -= weight(j);
    j++;
  }
  return j;
}



void SSA_direct::step(
    const NumericVector& state, 
    const NumericVector& transition_rates, 
    const NumericMatrix& nu,
    double* dtime, 
    NumericVector& dstate
) {
  int j = weighted_sample(transition_rates);
  
  for (int i = 0; i < dstate.size(); i++) {
    dstate(i) = nu(i, j);
  }
  
  *dtime = -log(runif(1, 0, 1)(0)) / sum(transition_rates); 
}



// // [[Rcpp::export]]
// List test_direct(NumericVector state, NumericVector transition_rates, NumericMatrix nu){
//   SSA_direct alg;
//   
//   double dtime;
//   NumericVector dstate(state.size());
//   // alg.step(state, transition_rates, nu, &dtime, &dstate);
//   alg.step(state, transition_rates, nu, &dtime, dstate);
//   
//   List ret; ret["dtime"] = dtime; ret["dstate"] = dstate;
//   return ret;
// }
