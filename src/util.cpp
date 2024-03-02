#include <numeric>
#include <math.h>
#include <algorithm>
#include <Rmath.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <RcppArmadillo.h>
// #include "dist_funs.h"

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat calc_distC(arma::mat X, arma::mat Y){
  int n_xrows = X.n_rows;
  int n_yrows = Y.n_rows;
  arma::mat D(n_xrows, n_yrows);
  for(int i = 0; i < n_xrows; i++){
    for(int j = 0; j < n_yrows; j++){
      D(i,j) = std::sqrt(arma::sum(pow(X.row(i) - Y.row(j), 2.0)));
    }
  }
  return(D);
}


// [[Rcpp::export]]
arma::mat calc_dist_selfC(arma::mat X){
  int n_xrows = X.n_rows;
  arma::mat D(n_xrows, n_xrows);
  for(int i = 0; i < n_xrows; i++){
    for(int j = i; j < n_xrows; j++){
      D(i,j) = std::sqrt(arma::sum(pow(X.row(i) - X.row(j), 2.0)));
      D(j,i) = D(i,j);
    }
  }
  return(D);
}

/*
 * Euclidean distance between two vectors
 */
// [[Rcpp::export]]
 double edistC(arma::vec x, arma::vec y){
   return(std::sqrt(arma::sum(pow(x - y, 2.0))));
 }

 // [[Rcpp::export]]
 double gauss_kernelC(arma::vec x, arma::vec y, double theta = 1.0){
   return( exp( -pow(edistC(x, y), 2.0) / (2.0 * theta))) ;
 }

 // [[Rcpp::export]]
 arma::mat gauss_kernel_gramC(arma::mat X, double theta = 1.0){
   int n_xrows = X.n_rows;
   arma::mat G(n_xrows, n_xrows);
   double tmp;
   for(int i = 0; i < n_xrows; i++){
     for(int j = i; j < n_xrows; j++){
       tmp = gauss_kernelC(arma::conv_to<arma::vec>::from(X.row(i)),
                           arma::conv_to<arma::vec>::from(X.row(j)),
                           theta);
       G(i,j) = tmp;
       G(j,i) = tmp;
     }
   }
   return(G);
 }

/*
* Transform correlated multivariate normal random vectors to iid normal random
* vectors with mean 0 and variance 1.
* Arguments:
* Z: p x n matrix of mean-zero multivariate normal samples (columns contain
*    independent samples)
* Sigma p x p covariance matrix
*
* Value:
* Matrix (p x n) of whitened Z values (iid, mean 0, variance 1)
*/
// [[Rcpp::export]]
arma::mat decorrelateC(arma::mat Z, arma::mat Sigma){
  arma::mat U(Sigma.n_rows, Sigma.n_cols);
  U = arma::chol(Sigma);
  int n = Z.n_cols;
  for(int i = 0; i < n; i++){
    Z.col(i) = arma::solve(U.t(), Z.col(i));
  }
  return(Z);
}

// [[Rcpp::export]]
arma::mat cholC(arma::mat Sigma){
  return(arma::chol(Sigma));
}

// [[Rcpp::export]]
arma::mat solveC(arma::mat Sigma){
  return(arma::inv_sympd(Sigma));
}

////[[Rcpp::export]]
//arma::vec pnormC(arma::vec x, arma::vec mu, arma::vec sd) // Phi(-∞, x) aka N(x)
//{
//  int n = x.n_elem;
//  arma::vec p(n);
//  for(int i = 0; i < n; i++){
//    p[i] = std::erfc(-((x[i] - mu[i])/sd[i])/std::sqrt(2))/2;
//  }
//  return(p);
//}

// [[Rcpp::export]]
Rcpp::NumericVector myqnorm(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma) {
  return Rcpp::qnorm((x - mu)/sigma);
}

/*
* Simulate from a truncated multivariate normal distribution using Gibbs sampling. This is a c++ port
* from the R package tmvtnorm::rtmvnorm.
* Arguments:
*   n: number of samples to draw
*   mean: vector of Normal means of non-truncated distribution (length = p)
*   H: precision matrix (inverse covariance matrix) dimension p x p
*   lower: vector of lower truncation bounds (length = p)
*   upper: vector of upper truncation bounds (length = p)
*   burn_in_samples: (integer) number of burn
* Value:
*   matrix (n x p dimensional) of truncated normal random samples
*
* Comments: If some components of lower and upper are the same (i.e. corresponding to a degenerate
*           distribution in that dimension, the simulated values will not be exactly equal
*           to that common value of lower and upper. They will be equal up to machine
*           precision, but there will be some noise.)
*
*/

// [[Rcpp::export]]
arma::mat rtmvtnormC (int n,
                    arma::vec mean,
                    arma::mat H,
                    arma::vec lower,
                    arma::vec upper,
                    arma::vec init,
                    int burn_in_samples = 0,
                    int thinning = 1,
                    bool initialize_x0 = false)
{
  if (thinning < 1){
    throw std::range_error("Thinning number must be positive");
  }
  int d = mean.n_elem;
  int S = burn_in_samples;
  if((lower.n_elem != d)|(upper.n_elem != d)){
    throw std::range_error("lower, upper, and mean vectors must all be of the same length");
  }
  if (S < 0){
    throw std::range_error("number of burn-in samples must be non-negative");
  }
  // Set start value for sample
  arma::vec x0(d);
  if(initialize_x0 == true){
                            // column vector
    Rcpp::LogicalVector linf = Rcpp::is_infinite(as<NumericVector>(wrap(lower)));
    Rcpp::LogicalVector uinf = Rcpp::is_infinite(as<NumericVector>(wrap(upper)));
    for(int i = 0; i < d; i++){
      if(!linf[i]){
        x0[i] = lower[i];
      }
      else{
        if(!uinf[i]){
          x0[i] = upper[i];
        }
        else{
          x0[i] = 0.0;
        }
      }
    }
  }
  else{
    x0 = init;
  }
  arma::mat X(n, d, arma::fill::randu);
  arma::vec U((S + n * thinning + 2) * d);
  U.fill(arma::fill::randu);
  int l = 1;
  arma::vec x = x0;
  arma::vec sd = sqrt(1/H.diag());
  arma::mat mu_i(1,1);
  arma::colvec h(d);
  arma::colvec h_m1(d-1);
  arma::colvec idx = arma::linspace<arma::vec>(0, d-1, d);
  arma::uvec minus_i(d-1);
  double Fa;
  double Fb;
  double za;
  double zb;
  double qn;

  for(int j = -S; j < n * thinning; j++){
    for(int i = 0; i < d; i ++){
      minus_i= find(idx != i);
      h = H.col(i);
      h_m1 = h.elem(minus_i);
      mu_i = mean[i] - (1/H(i,i)) *  (h_m1.t() * (x.elem(minus_i) - mean.elem(minus_i)));
      za = Rcpp::as<double>(wrap((lower[i]-mu_i)/sd[i]));
      zb = Rcpp::as<double>(wrap((upper[i]-mu_i)/sd[i]));
      Fa = Rcpp::as<double>(wrap(Rcpp::pnorm(Rcpp::as<NumericVector>(wrap(za)))));
      Fb = Rcpp::as<double>(wrap(Rcpp::pnorm(Rcpp::as<NumericVector>(wrap(zb)))));
      qn = U[l] * (Fb - Fa) + Fa;
      x[i] = Rcpp::as<double>(wrap(mu_i + sd[i] * Rcpp::as<double>(wrap(Rcpp::qnorm(Rcpp::as<NumericVector>(wrap(qn)))))));
      l++;

    }
    if(j > -1){
      if(thinning == 1){
        X.row(j) = x.t();
      }
      else if( j % thinning == 0){
        X.row(floor(j/thinning)) = x.t();
      }
    }
  }
  return(X);
}


/*
 * Matrix version of rtmvtnorm. That can generate independent, but not
 * identically distributed multivariate truncated normals. In the notation below
 * n denotes the sample number where, e.g., different samples can come from
 * multivariate truncated normal distributions with different means, but with a
 * common covariance matrix. The d
 * Arguments:
 *   mean: Matrix (n x p) of Normal means of non-truncated distribution (length = p)
 *   H: precision matrix (inverse covariance matrix) dimension p x p. This is
 *      common for the distributions of all n samples
 *   lower: matrix (n x p) of lower truncation bounds
 *   upper: matrix (n x p) of upper truncation bounds
 *   burn_in_samples: (integer) number of burn
 * Value:
 *   matrix (n x p dimensional) of truncated normal random samples
 *
 * Comments: If some components of lower and upper are the same (i.e. corresponding to a degenerate
 *           distribution in that dimension, the simulated values will not be exactly equal
 *           to that common value of lower and upper. They will be equal up to machine
 *           precision, but there will be some noise.)
 *
 */

// // [[Rcpp::export]]
// arma::mat MrtmvtnormC ( // int n,
//                     arma::mat mean,
//                     arma::mat H,
//                     arma::mat lower,
//                     arma::mat upper,
//                     int burn_in_samples = 0,
//                     int thinning = 1){
//
//   int nrep = mean.n_cols;
//   int nloc = mean.n_rows;
//   arma::mat z(nloc, nrep);
//   for(int i = 0; i < nrep; i++){
//     z.col(i) = rtmvtnormC (1,
//                           mean.col(i),
//                           H,
//                           lower.col(i),
//                           upper.col(i),
//                           burn_in_samples,
//                           thinning).t();
//    }
//   return(z);
// }

