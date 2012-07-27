#include <Rcpp.h>
#include <RcppArmadillo.h>

extern "C" SEXP bbFixefGroup(SEXP groupIndex, SEXP X_, SEXP Vinv_, 
                             SEXP XVXinv_, SEXP P_, SEXP e_) {

  Rcpp::List index(groupIndex);
  arma::mat X = Rcpp::as<arma::mat>(X_);
//  arma::mat Zt = Rcpp::as<arma::mat>(Zt_);
  arma::mat XVXinv = Rcpp::as<arma::mat>(XVXinv_);
  arma::mat Vinv = Rcpp::as<arma::mat>(Vinv_);
//  arma::mat D = Rcpp::as<arma::mat>(D_);
  arma::mat P = Rcpp::as<arma::mat>(P_);
  arma::colvec e = Rcpp::as<arma::colvec>(e_);
  
  int n = index.size();
  arma::mat Xt;
  arma::mat Vi, Pii, Pi, beta_cdd;
  arma::colvec ei;
  Rcpp::List reslist( n );
  int ii;

  for(ii=0; ii<n; ii++){
     Rcpp::IntegerVector tindex = index[ii];
     arma::uvec tind = Rcpp::as<arma::uvec>(tindex);
     Vi  = Vinv.cols(tind);
     Pii = P(tind, tind);
//     Pi  = P.cols(tind);
     ei  = e.elem(tind);
     
    beta_cdd  = XVXinv * trans(X) * Vi * inv(Pii) * ei; // beta - beta_(i)
          
//     reslist[ii] = Rcpp::List::create(Rcpp::Named("beta_cdd") = beta_cdd);
     reslist[ii] = beta_cdd;
  }
  return reslist;
}