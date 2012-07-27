#include <Rcpp.h>
#include <RcppArmadillo.h>

extern "C" SEXP bbRanef(SEXP Zt_, SEXP D_, SEXP P_, SEXP e_) {
  arma::mat Zt = Rcpp::as<arma::mat>(Zt_);
  arma::mat D  = Rcpp::as<arma::mat>(D_);
  arma::mat P  = Rcpp::as<arma::mat>(P_);
  arma::colvec e = Rcpp::as<arma::colvec>(e_);
  
  int n = Zt.n_cols;
  arma::colvec b_cdd, ei;
  Rcpp::List reslist( n );
  int ii;
  
  for(ii=0; ii<n; ii++){
    ei = arma::zeros( n ); 
    ei(ii) = e(ii);
    b_cdd = D * Zt * P * ei;
    reslist[ii] = b_cdd;
  }
  
  return reslist;
}