#ifndef DYNGEN_FASTGSSA_H
#define DYNGEN_FASTGSSA_H

#include <Rcpp.h>
#include "ssa.hpp"
#include "ssa_em.hpp"
#include "ssa_direct.hpp"

using namespace Rcpp;

// RCPP_EXPOSED_CLASS(SSA)
// RCPP_EXPOSED_CLASS(SSA_EM)
// RCPP_EXPOSED_CLASS(SSA_direct)
// 
// RCPP_MODULE(fastgssa){
//   Rcpp::class_<SSA>("SSA")
//   .method("simulate", &SSA::simulate)
//   ;
//   Rcpp::class_<SSA_EM>("SSA_EM")
//     .derives<SSA>("SSA")
//     .constructor<double, double>()
//     .field("h", &SSA_EM::h)
//     .field("noise_strength", &SSA_EM::noise_strength)
//   ;
//   Rcpp::class_<SSA_direct>("SSA_direct")
//     .derives<SSA>("SSA")
//     .constructor()
//   ;
// }

#endif