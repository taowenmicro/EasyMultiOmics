#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @title Schur Product
 //' @description Element-wise multiplication of two matrices (Hadamard product)
 //' @param a Matrix A
 //' @param b Matrix B
 //' @return Matrix
 //' @keywords internal
 // [[Rcpp::export]]
 arma::mat schur(arma::mat& a, arma::mat& b) {
   return(a % b);
 }
