#ifndef DYNGEN_SSA_DIRECT_H
#define DYNGEN_SSA_DIRECT_H

#include <Rcpp.h>
#include "ssa.h"

using namespace Rcpp;

class SSA_direct : public SSA {
public:
  SSA_direct() : SSA() {} 
  void step(
      const NumericVector& state, 
      const NumericVector& transition_rates, 
      const NumericMatrix& nu,
      double* dtime, 
      NumericVector& dstate
  );
} ;

#endif
