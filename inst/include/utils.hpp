#ifndef DYNGEN_UTILS_H
#define DYNGEN_UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

List resize( const List& x, int n ){
  int oldsize = x.size() ;
  List y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}

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


#endif