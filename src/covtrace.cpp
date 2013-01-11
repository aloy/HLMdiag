#include <RcppArmadillo.h>

extern "C" SEXP covtraceCalc(SEXP index, SEXP X_, SEXP P_,
							  SEXP Vinv_, SEXP XVXinv_) {
  Rcpp::List grpInd(index);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat P = Rcpp::as<arma::mat>(P_);
  arma::mat Vinv = Rcpp::as<arma::mat>(Vinv_);
  arma::mat XVXinv = Rcpp::as<arma::mat>(XVXinv_);

  arma::mat Xt = trans(X);
  arma::mat XVX = inv(XVXinv);

  int n = grpInd.size(), p = X.n_cols;
  int ii;
  arma::mat Pa, XNaXinv;
  Rcpp::NumericVector covtrace ( n );

  for(ii=0; ii<n; ii++){
    arma::uvec ind = grpInd[ii];

    Pa = P.submat(ind, ind);
    XNaXinv = XVXinv + XVXinv * Xt * Vinv.cols(ind) * inv(Pa) * Vinv.rows(ind) * X * XVXinv;

    covtrace[ii] = trace(XVX * XNaXinv) - p;
  }
  return wrap( abs( covtrace ) );
}