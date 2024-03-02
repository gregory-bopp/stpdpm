#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector mypnorm(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma) {
  return Rcpp::pnorm((x - mu)/sigma);
}
