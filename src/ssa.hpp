#ifndef DYNGEN_SSA_H
#define DYNGEN_SSA_H

#include <Rcpp.h>
#include "utils.hpp"

using namespace Rcpp;

class SSA {
public:
  SSA(std::string name_) : name(name_) {}
  std::string name;
  
  virtual ~SSA() {}

  virtual void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const NumericMatrix& nu,
      NumericVector& dtime,
      NumericVector& dstate
  ) {
    stop("step() should have been overridden but wasn't!!!");
  }
} ;


#endif