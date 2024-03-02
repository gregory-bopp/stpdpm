# Extended GPD
# pegpd <- function(y, sigma, xi, kappa) {
#   return((1 - (pmax(1 + (
#     xi / sigma
#   ) * y, 0)) ^ (-1 / xi)) ^ kappa)
# }
#
# qegpd <- function(p, sigma, xi, kappa) {
#   return((sigma / xi) * (((1 - p ^ {
#     1 / kappa
#   }) ^ {
#     -xi
#   }) - 1))
# }
#
# degpd <- function(y, sigma, xi, kappa) {
#   return(((1 / sigma) * pmax(1 + (xi / sigma) * y, 0) ^ (-1 - 1 / xi)) *
#            (kappa * (1 - pmax(1 + xi / sigma * y, 0) ^ (-1 / xi)) ^ (kappa - 1)))
# }



#' Map from real line to positive reals that preserves right tail
#'
#' @param x input from real line
#' @param lambda skewness parameter (power-transform on positive real line
#' scale)
#'
#' @return mapped input on positive reals
#' @export
#'
#' @examples
#' x <- rnorm(1000)
#' y <- real_to_pos(x, lambda = 1)
#' par(mfrow = c(1,2))
#' hist(x, breaks = 10)
#' hist(y, breaks = 10)
#' curve(real_to_pos(x, lambda = 1/8), from = -1, to = 1, n = 10000)
real_to_pos <- function(x, lambda = 1) {
  log(1 + exp(x))^(1/lambda)
}


#' Map from real line to positive reals that preserves tails
#'
#' @param x input from the positive reals
#' @param lambda skewness parameter
#' @export
#' @examples
#' y <- rgamma(1000, shape = 1, scale = 30)
#' ps <- rep(NA, 100)
#' lambdas <- seq(0.1, 10, l = 100)
#' for(i in 1:length(lambdas)){
#'   x <- pos_to_real(y, lambda= lambdas[i])
#'   ps[i] <- shapiro.test(x)$p.value
#' }
#' lambda <- lambdas[which.max(ps)]
#' x <- pos_to_real(y, lambda= lambda)
#' par(mfrow = c(1,2))
#' hist(y, breaks = 10)
#' hist(x, breaks = 10)
#' curve(pos_to_real(x, lambda = 1/8), from = 0, to = 10, n = 10000)
pos_to_real <- function(x, lambda = 1) {
  log(exp(x ^ lambda) - 1)
}
