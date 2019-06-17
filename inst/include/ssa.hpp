#ifndef DYNGEN_SSA_H
#define DYNGEN_SSA_H

#include <Rcpp.h>
#include "utils.hpp"

using namespace Rcpp;

class SSA {
public:
  SSA() {}
  // SSA(TR_FUN fun) : calculate_transition_rates(fun) {}
  
  // TR_FUN calculate_transition_rates;
  // SSA(std::string name_) : name(name_) {}
  // std::string name;
  virtual ~SSA() {}

  virtual void step(
      const NumericVector& state,
      const NumericVector& transition_rates,
      const NumericMatrix& nu,
      NumericVector& dtime,
      NumericVector& dstate
  )
  // ) = 0;
  {
    Rcout << "step() should have been overridden but wasn't!!!" << std::endl;
    // this function should be overridden
  }
  
  // void calculate_transition_rates(
  //     const NumericVector& state,
  //     const NumericVector& params,
  //     const double time,
  //     NumericVector& transition_rates
  // ) {
  //   Rcout << "calculate_transition_rates() should have been overridden but wasn't!!!" << std::endl;
  // };
  

  // virtual void calculate_transition_rates(
  //     const NumericVector& state,
  //     const NumericVector& params,
  //     const double time,
  //     NumericVector& transition_rates
  // ) = 0;
  // {
  //   Rcout << "calculate_transition_rates() should have been overridden but wasn't!!!" << std::endl;
  //   // this function should be overridden
  // }

} ;


#endif