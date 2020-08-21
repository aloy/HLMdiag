#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

extern "C" SEXP mdffitsSubset(SEXP index, SEXP X_, SEXP P_,
							   SEXP Vinv_, SEXP XVXinv_, SEXP e_) {
  Rcpp::List grpInd(index);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat P = Rcpp::as<arma::mat>(P_);
  arma::mat Vinv = Rcpp::as<arma::mat>(Vinv_);
  arma::mat XVXinv = Rcpp::as<arma::mat>(XVXinv_);
  arma::colvec e = Rcpp::as<arma::colvec>(e_);

  arma::mat Xt = trans(X);

  int n = grpInd.size(), p = X.n_cols;
  int ii;
  arma::mat Pa, Na, XNaXinv, XNaX;
  arma::colvec cdd;
  Rcpp::NumericVector mdffits ( n );
  Rcpp::List beta_cdd ( n );

  for(ii=0; ii<n; ii++){
    arma::uvec ind = grpInd[ii];

    Pa = P.submat(ind, ind);
/*
    Na = Vinv - Vinv.cols(ind) * inv( Vinv.submat(ind, ind) ) * Vinv.rows(ind);
    XNaX = Xt * Na * X;
*/
	
	XNaXinv = XVXinv + XVXinv * Xt * Vinv.cols(ind) * inv(Pa) * Vinv.rows(ind) * X * XVXinv;
	XNaX = inv(XNaXinv);

    cdd = XVXinv * Xt * Vinv.cols(ind) * inv(Pa) * Vinv.rows(ind) * e;
    mdffits[ii] = arma::as_scalar( trans(cdd) * XNaX * cdd ) / p;
    beta_cdd[ii] = cdd;

  }

  return Rcpp::List::create( Rcpp::Named("mdffits") = mdffits, 
                             Rcpp::Named("beta_cdd") = beta_cdd);
}