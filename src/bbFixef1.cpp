#include <Rcpp.h>
#include <RcppArmadillo.h>

extern "C" SEXP bbFixef1(SEXP y_, SEXP X_, SEXP Vinv_, SEXP XVXinv_, SEXP beta_) {
  arma::colvec Y = Rcpp::as<arma::colvec>(y_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Vinv = Rcpp::as<arma::mat>(Vinv_);
  arma::mat XVXinv = Rcpp::as<arma::mat>(XVXinv_);
  arma::colvec beta = Rcpp::as<arma::colvec>(beta_);

  int n = Vinv.n_rows;
  double sbb, hbb, ybb;
  arma::mat Xt;
  arma::colvec Vi, xbb, beta_cdd;
  Rcpp::List reslist( n );
  int ii;

  for(ii=0; ii<n; ii++){
     Vi = Vinv.col(ii);
         
     sbb = 1 / Vinv(ii,ii);
     xbb = sbb * trans(X) * Vi;
     hbb = arma::as_scalar(trans(xbb) * XVXinv * xbb);
     ybb = arma::as_scalar(sbb * trans(Y) * Vi);
     
     beta_cdd = beta - (1 / (sbb - hbb)) * XVXinv * xbb * (ybb - trans(xbb) * beta); // beta_(i)
          
     reslist[ii] = Rcpp::List::create(Rcpp::Named("sbb") = sbb, Rcpp::Named("xbb") = xbb,
                                      Rcpp::Named("hbb") = hbb, Rcpp::Named("ybb") = ybb,
                                      Rcpp::Named("beta_cdd") = beta_cdd);
  }
  return reslist;
}