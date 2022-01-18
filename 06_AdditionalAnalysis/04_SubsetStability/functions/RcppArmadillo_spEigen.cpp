#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List sp_eigen_cpp(arma::sp_mat M,
                        int k) {
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  /* arma::eig_gen(eigval, eigvec, M, "balance"); */
  arma::eigs_gen(eigval, eigvec, M, k, "lm");
  
  return Rcpp::List::create(eigval, eigvec);
}
