#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

SEXP cooksdSubset(SEXP index, SEXP X_, SEXP P_, 
                              SEXP Vinv_, SEXP XVXinv_, SEXP e_) {
  Rcpp::List grpInd(index);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat P = Rcpp::as<arma::mat>(P_);
  arma::mat Vinv = Rcpp::as<arma::mat>(Vinv_);
  arma::mat XVXinv = Rcpp::as<arma::mat>(XVXinv_);
  arma::colvec e = Rcpp::as<arma::colvec>(e_);

  arma::mat Xt = trans(X), XVX = inv(XVXinv);

  int n = grpInd.size(), p = X.n_cols;
  int ii;
  arma::mat Pa;
  arma::colvec cdd;
  Rcpp::NumericVector cooksd ( n );
  Rcpp::List beta_cdd ( n );

  for(ii=0; ii<n; ii++){
    arma::uvec ind = grpInd[ii];

    Pa = P.submat(ind, ind);

    cdd = XVXinv * Xt * Vinv.cols(ind) * inv(Pa) * Vinv.rows(ind) * e;
    cooksd[ii] = arma::as_scalar( trans(cdd) * XVX * cdd ) / p;
    beta_cdd[ii] = cdd;

  }

  return Rcpp::List::create( Rcpp::Named("cooksd") = cooksd, 
                             Rcpp::Named("beta_cdd") = beta_cdd);
}