#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

extern "C" SEXP cxxmatsub(SEXP BB, SEXP CC) {
  arma::mat B = Rcpp::as<arma::mat>(BB);
  arma::mat C = Rcpp::as<arma::mat>(CC);
  return Rcpp::wrap(B - C);
}