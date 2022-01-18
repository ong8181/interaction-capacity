#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List eigen_cpp(arma::fmat M) {
  arma::cx_fvec eigval;
  arma::cx_fmat eigvec;
  /* arma::eig_gen(eigval, eigvec, M, "balance"); */
  arma::eig_gen(eigval, eigvec, M);
  
  return Rcpp::List::create(eigval, eigvec);
}
