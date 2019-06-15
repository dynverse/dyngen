#include <Rcpp.h>
#include "ssa_em.h"

using namespace Rcpp;

void SSA_EM::step(
    const NumericVector& state, 
    const NumericVector& transition_rates, 
    const NumericMatrix& nu,
    double* dtime, 
    NumericVector& dstate
) {
  *dtime = h;
  
  for (int i = 0; i < dstate.size(); i++) {
    double out = 0.0;
    for (int j = 0; j < transition_rates.size(); j++) {
      out += nu(i, j) * transition_rates(j) * h;
    }
    out += sqrt(abs(state(i))) * noise_strength * rnorm(1, 0.0, h)(0);
    dstate(i) = out;
  }
}

// // [[Rcpp::export]]
// List test_em(double h, double noise_strength, NumericVector state, NumericVector transition_rates, NumericMatrix nu){
//   SSA_EM alg(h, noise_strength);
//   
//   double dtime;
//   NumericVector dstate(state.size());
//   alg.step(state, transition_rates, nu, &dtime, dstate);
//   
//   List ret; ret["dtime"] = dtime; ret["dstate"] = dstate;
//   return ret;
// }

RCPP_EXPOSED_CLASS(SSA_EM)
RCPP_MODULE(SSA_EM){
  Rcpp::class_<SSA>("SSA")
  .method("simulate", &SSA::simulate)
  ;
  Rcpp::class_<SSA_EM>("SSA_EM")
    .derives<SSA>("SSA")
    .constructor<double, double>()
    .field("h", &SSA_EM::h) 
    .field("noise_strength", &SSA_EM::noise_strength) 
  ;
}
