library(testthat)
library(stpdpm)

context("Testing utility functions")

test_edist <- function(){
  x <- exp(rnorm(100, sd = 10))
  y <- exp(rnorm(100, sd = 10))
  e1 <- edist(x,y)
  e2 <- edistC(x,y)
  expect_equal(e2, e1)
}

test_calc_dists <- function(){
  x <- matrix(rnorm(100^2, sd= 1000), nrow = 20, ncol = 20)
  y <- matrix(rnorm(100^2, sd = 1000), nrow = 20, ncol = 20)
  d1 <- calc_distR(x,y)
  d2 <- calc_dist(x,y)
  expect_equal(d1, d2)
}

test_gauss_kernel <- function(){
  x <- rnorm(1e2, sd = 1)
  y <- rnorm(1e2, sd = 1)
  e1 <- gauss_kernel(calc_dist(matrix(x, ncol =  1), matrix(x, ncol = 1)), theta = 1)
  e2 <- gauss_kernel_gram(matrix(x,ncol = 1), theta = 1)
  expect_equal(e2, e1)
}


test_gauss_kernel_gram <- function(){
  X <- matrix(rnorm(200*10), nrow = 1000, ncol = 10)
  G <- gauss_kernel_gramC(X, 1)
  GR <- gauss_kernel_gram(X, 1)
  expect_equal(G, GR)
}

# Make sure inverse gamma agrees with dgamma
test_dinvgamma <- function(){
  a <- rexp(1)
  b <- rexp(1)
  y<- 1/rgamma(100000, shape = a, rate = b)
  val = dinvgamma(y, shape = a, scale = b, log = T)
  expected = dgamma(1/y, shape = a, scale = 1/b, log = T) - 2*log(y)
  expect_equal(val,expected)
}

test_dmvnp <- function(){
  library(mvtnorm)
  n <- 100
  p <- 10
  A <- matrix(runif(p^2), nrow = p, ncol= p)
  A <- (t(A) + A)/2 + diag(p)*p
  iA <- solve(A)
  mu <- rnorm(p)
  x <- rmvnorm(n, mean = mu, sigma = A)
  expected = dmvnorm(x, mu, A, log = T)
  value = dmvnp(x, mu, iA, log = T)
  expect_equal(value, expected)
}

# Test that the theoretical mean is close to the simulated one.
# This is a weak test that may fail even if the code is okay.
test_rtmvtnormC <- function(){
  set.seed(2)
  n <- 10000
  p <- 10
  A <- matrix(runif(p^2), nrow = p, ncol= p)
  A <- (t(A) + A)/2 + diag(p)*p
  iA <- solve(A)
  mu <- rnorm(p)
  y = rtmvtnormC(n,
                 mean = mu,
                 H = iA,
                 lower = rep(0,p),
                 upper = rep(Inf, p),
                 burn_in_samples = 10,
                 thinning = 1,
                 init = 0,
                 initialize_x0 = T)
  alpha <- (0 - mu)/sqrt(diag(A))
  expected <- mu + dnorm(alpha)*sqrt(diag(A))/(1-pnorm(alpha))
  expect_true(abs(mean(apply(y,2, mean) - expected)) < 0.5)
}

test_is_scalar <- function(){
  expect_true(is_scalar(3))
}


test_that("Testing distance and kernel functions", {
  test_edist()
  test_calc_dists()
  test_gauss_kernel()
  test_gauss_kernel_gram()
  test_dinvgamma()
  test_dmvnp()
  test_rtmvtnormC()
  test_is_scalar()
})
