#ifndef DYNGEN_SSA_EM_H
#define DYNGEN_SSA_EM_H

#include <Rcpp.h>
#include "ssa.hpp"

using namespace Rcpp;

class SSA_EM : public SSA {
public:
   SSA_EM(double h_, double noise_strength_) : SSA("EM"), h(h_), noise_strength(noise_strength_) {}
  // SSA_EM(double h_, double noise_strength_) : SSA(), h(h_), noise_strength(noise_strength_) {}
  // SSA_EM(TR_FUN fun, double h_, double noise_strength_) : SSA(fun), h(h_), noise_strength(noise_strength_) {}

  double h ;
  double noise_strength ;

  void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const NumericMatrix& nu,
      NumericVector& dtime,
      NumericVector& dstate
  ) {
    dtime(0) = h;

    for (int i = 0; i < dstate.size(); i++) {
      double out = 0.0;
      for (int j = 0; j < transition_rates.size(); j++) {
        out += nu(i, j) * transition_rates(j) * h;
      }
      out += sqrt(abs(state(i))) * noise_strength * rnorm(1, 0.0, h)(0);
      dstate(i) = out;
    }
  }
} ;

//' @export
// [[Rcpp::export]]
SEXP ssa_em(double h = .01, double noise_strength = 2.0) {
  SSA_EM *ssa = new SSA_EM(h, noise_strength);
  XPtr<SSA_EM> ptr(ssa);
  return ptr;
}

#endif
