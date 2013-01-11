#include <RcppArmadillo.h>

extern "C" SEXP cooksdObs(SEXP y_, SEXP X_, SEXP Vinv_, 
                           SEXP XVXinv_, SEXP beta_) {
  arma::colvec Y = Rcpp::as<arma::colvec>(y_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Vinv = Rcpp::as<arma::mat>(Vinv_);
  arma::mat XVXinv = Rcpp::as<arma::mat>(XVXinv_);
  arma::colvec beta = Rcpp::as<arma::colvec>(beta_);

  int p = X.n_cols;
  int n = Vinv.n_rows;
  double sbb, hbb, ybb;
  arma::mat Xt = trans(X), XVX = inv(XVXinv);
  arma::colvec Vi, xbb, cdd;
  Rcpp::NumericVector cooksd ( n );
  Rcpp::List beta_cdd( n );
  int ii;

  for(ii=0; ii<n; ii++){
     Vi = Vinv.col(ii);

     sbb = 1 / Vinv(ii,ii);
     xbb = sbb * trans(X) * Vi;
     hbb = arma::as_scalar(trans(xbb) * XVXinv * xbb);
     ybb = arma::as_scalar(sbb * trans(Y) * Vi);

     cdd = (1 / (sbb - hbb)) * XVXinv * xbb * (ybb - trans(xbb) * beta);
     cooksd[ii] = arma::as_scalar( trans(cdd) * XVX * cdd) / p;
     beta_cdd[ii] = cdd;
  }
  return Rcpp::List::create( Rcpp::Named("cooksd") = cooksd, 
                             Rcpp::Named("beta_cdd") = beta_cdd);
}