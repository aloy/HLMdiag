#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

extern "C" SEXP covratioCalc(SEXP index, SEXP X_, SEXP P_,
							  SEXP Vinv_, SEXP XVXinv_) {
	Rcpp::List grpInd(index);
	arma::mat X = Rcpp::as<arma::mat>(X_);
	arma::mat P = Rcpp::as<arma::mat>(P_);
	arma::mat Vinv = Rcpp::as<arma::mat>(Vinv_);
	arma::mat XVXinv = Rcpp::as<arma::mat>(XVXinv_);
	
	arma::mat Xt = trans(X);
	
	int n = grpInd.size();
	int ii;
	arma::mat Pa, XNaXinv;
	Rcpp::NumericVector covratio ( n );
	
	for(ii=0; ii<n; ii++){
		arma::uvec ind = grpInd[ii];
		
		Pa = P.submat(ind, ind);
		XNaXinv = XVXinv + XVXinv * Xt * Vinv.cols(ind) * inv(Pa) * Vinv.rows(ind) * X * XVXinv;
		
		covratio[ii] =  det(XNaXinv) / det(XVXinv);
	}	
	return wrap( covratio );
}