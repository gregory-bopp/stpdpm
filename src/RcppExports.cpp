// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// mypnorm
Rcpp::NumericVector mypnorm(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma);
RcppExport SEXP _stpdpm_mypnorm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mypnorm(x, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// calc_distC
arma::mat calc_distC(arma::mat X, arma::mat Y);
RcppExport SEXP _stpdpm_calc_distC(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_distC(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// calc_dist_selfC
arma::mat calc_dist_selfC(arma::mat X);
RcppExport SEXP _stpdpm_calc_dist_selfC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_dist_selfC(X));
    return rcpp_result_gen;
END_RCPP
}
// edistC
double edistC(arma::vec x, arma::vec y);
RcppExport SEXP _stpdpm_edistC(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(edistC(x, y));
    return rcpp_result_gen;
END_RCPP
}
// gauss_kernelC
double gauss_kernelC(arma::vec x, arma::vec y, double theta);
RcppExport SEXP _stpdpm_gauss_kernelC(SEXP xSEXP, SEXP ySEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_kernelC(x, y, theta));
    return rcpp_result_gen;
END_RCPP
}
// gauss_kernel_gramC
arma::mat gauss_kernel_gramC(arma::mat X, double theta);
RcppExport SEXP _stpdpm_gauss_kernel_gramC(SEXP XSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_kernel_gramC(X, theta));
    return rcpp_result_gen;
END_RCPP
}
// decorrelateC
arma::mat decorrelateC(arma::mat Z, arma::mat Sigma);
RcppExport SEXP _stpdpm_decorrelateC(SEXP ZSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(decorrelateC(Z, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// cholC
arma::mat cholC(arma::mat Sigma);
RcppExport SEXP _stpdpm_cholC(SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(cholC(Sigma));
    return rcpp_result_gen;
END_RCPP
}
// solveC
arma::mat solveC(arma::mat Sigma);
RcppExport SEXP _stpdpm_solveC(SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(solveC(Sigma));
    return rcpp_result_gen;
END_RCPP
}
// myqnorm
Rcpp::NumericVector myqnorm(Rcpp::NumericVector x, Rcpp::NumericVector mu, Rcpp::NumericVector sigma);
RcppExport SEXP _stpdpm_myqnorm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(myqnorm(x, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// rtmvtnormC
arma::mat rtmvtnormC(int n, arma::vec mean, arma::mat H, arma::vec lower, arma::vec upper, arma::vec init, int burn_in_samples, int thinning, bool initialize_x0);
RcppExport SEXP _stpdpm_rtmvtnormC(SEXP nSEXP, SEXP meanSEXP, SEXP HSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP initSEXP, SEXP burn_in_samplesSEXP, SEXP thinningSEXP, SEXP initialize_x0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type burn_in_samples(burn_in_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type thinning(thinningSEXP);
    Rcpp::traits::input_parameter< bool >::type initialize_x0(initialize_x0SEXP);
    rcpp_result_gen = Rcpp::wrap(rtmvtnormC(n, mean, H, lower, upper, init, burn_in_samples, thinning, initialize_x0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stpdpm_mypnorm", (DL_FUNC) &_stpdpm_mypnorm, 3},
    {"_stpdpm_calc_distC", (DL_FUNC) &_stpdpm_calc_distC, 2},
    {"_stpdpm_calc_dist_selfC", (DL_FUNC) &_stpdpm_calc_dist_selfC, 1},
    {"_stpdpm_edistC", (DL_FUNC) &_stpdpm_edistC, 2},
    {"_stpdpm_gauss_kernelC", (DL_FUNC) &_stpdpm_gauss_kernelC, 3},
    {"_stpdpm_gauss_kernel_gramC", (DL_FUNC) &_stpdpm_gauss_kernel_gramC, 2},
    {"_stpdpm_decorrelateC", (DL_FUNC) &_stpdpm_decorrelateC, 2},
    {"_stpdpm_cholC", (DL_FUNC) &_stpdpm_cholC, 1},
    {"_stpdpm_solveC", (DL_FUNC) &_stpdpm_solveC, 1},
    {"_stpdpm_myqnorm", (DL_FUNC) &_stpdpm_myqnorm, 3},
    {"_stpdpm_rtmvtnormC", (DL_FUNC) &_stpdpm_rtmvtnormC, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_stpdpm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
