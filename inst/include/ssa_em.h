#ifndef DYNGEN_SSA_EM_H
#define DYNGEN_SSA_EM_H

#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_EM : public SSA {
public:
  SSA_EM(double h_, double noise_strength_) : SSA(), h(h_), noise_strength(noise_strength_) {} 
  double h ;
  double noise_strength ;
  void step(
      const NumericVector& state, 
      const NumericVector& transition_rates, 
      const NumericMatrix& nu,
      double* dtime, 
      NumericVector& dstate
  );
} ;

#endif
