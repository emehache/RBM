#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double distPaQ_cpp(NumericMatrix P,NumericMatrix Q){
  NumericVector output(P.nrow());
  for(int i=0;i<P.nrow();i++){
    output[i]=min( pow(P(i,0)-Q(_,0),2) + pow(P(i,1)-Q(_,1),2));
      }
  return sqrt(max(output)) ;
}
