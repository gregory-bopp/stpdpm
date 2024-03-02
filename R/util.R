#' @useDynLib stpdpm
dummy <- function(x){return(x)}


#' Check that object is a scalar
#'
#' @param x input object to test
#'
#' @return (logical) is x an atomic type of length 1
#' @export
#' @examples
#' is_scalar(3)  # TRUE
#' is_scalar(c(3,4)) # FALSE
is_scalar <- function(x){ is.atomic(x) && length(x) == 1L}

#' Repeat scalars and leave vectors alone
#' @param x scalar or vector
#' @param times number of times to repeat x if it is a scalar
#'
#' @export
#' @examples
#' expand_one_scalar(3,2)
#' expand_one_scalar(c(1,2), 5)
expand_one_scalar <- function(x, times) {
  if (is_scalar(x)) {
    return(rep(x, times = times))
  }
  else{
    return(x)
  }
}

#' Repeat scalars and leave vectors alone
#' @description This function modifies variables in its parent environment
#' but does not return anything. It replaces scalars by vectors of that scalar
#' repeated \code{reps} times. If the variable is already a vector it is left
#' unchanged.
#' @param varnames a vector of variable name strings
#' @param reps number of times the variables names in \code{varnames} that
#' are scalars should be repeated
#'
#' @export
#' @examples
#' a <- 1
#' b <- c(1,2,3)
#' expand_scalars(c("a","b"), 3)  # should repeat a 3 times and leave b alone
#' print(a); print(b)
expand_scalars <- function(varnames, reps) {
  for (i in 1:length(varnames)) {
    assign(varnames[i],
           do.call(expand_one_scalar, args = list(
             x = eval(parse(text = varnames[i]),
                      envir = parent.frame()),
             times = reps
           )),
           pos = parent.frame())
  }
  return(NULL)
}



#' Euclidean distance
#'
#' @param x vector 1
#' @param y vector 2
#'
#' @return Euclidean distance between x and y
#' @examples
#' x <- c(0,0)
#' y <- c(1,1)
#' edist(x, y)
edist <- function(x, y) {
  return(sqrt(sum((x - y) ^ 2)))
}


#' Gaussian kernel function
#'
#' @param dists matrix of distances
#' @param theta kernel scale parameter (i.e. variance in gauss density)
#'
#' @export
gauss_kernel <- function(dists, theta = 1) {
  return(exp(-(dists ^ 2) / (2 * theta)))
}


#' Multivariate Inverse logit function
#'
#' @param alpha (vector) of reals with first element = 0
#' @description the first element is treated as the reference category
#' @return vector of weights (0 < w < 1) that sum to 1
#' @export
#'
#' @examples
#' a <- c(0, rnorm(5))
#' ilogit(a)
ilogit <- function(alpha) {
  return(exp(alpha) / sum(exp(alpha)))
}

#' Multivariate logit function
#' @description the first element is treated as the reference category
#' @param p vector of weights each between 0 and 1
#'
#' @return vector of real numbers (first element is 0)
#' @export
#'
#' @examples
#' a <- c(0, rnorm(5))
#' max(logit(ilogit(a)) - a)
logit <- function(p) {
  c(0, log(p[-1] / (p[1])))
}



#' Whitened Multivariate normal vectors
#'
#' @param Y p x n matrix of observations (columns contain independent
#' replicates)
#' @param S p x p covariance matrix
#'
#' @return i.i.d. version of Y on N(0,1)
#' @export
Rdecorrelate <- function(Y, S) {
  U <- chol(S)
  forwardsolve(t(U), Y)
}


#' @export
rtmvtnormR <-
  function (n,
            mean = rep(0, nrow(H)),
            H = diag(length(mean)),
            lower = rep(-Inf, length = length(mean)),
            upper = rep(Inf, length = length(mean)),
            burn_in_samples = 0,
            start_value = NULL,
            thinning = 1)
  {
    if (thinning < 1 || !is.numeric(thinning) || length(thinning) >
        1) {
      stop("thinning must be a integer scalar > 0")
    }
    d <- length(mean)
    S <- burn_in_samples
    if (!is.null(S)) {
      if (S < 0)
        stop("number of burn-in samples must be non-negative")
    }
    if (!is.null(start_value)) {
      if (length(mean) != length(start_value))
        stop("mean and start value have non-conforming size")
      if (any(start_value < lower || start_value > upper))
        stop("start value is not inside support region")
      x0 <- start_value
    }
    else {
      x0 <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper),
                                                   upper, 0))
    }
    if (d == 1) {
      X <- rtnorm.gibbs(
        n,
        mu = mean[1],
        sigma = 1 / H[1, 1],
        a = lower[1],
        b = upper[1]
      )
      return(X)
    }
    X <- matrix(NA, n, d)
    U <- runif((S + n * thinning) * d)
    l <- 1
    sd <- sqrt(1 / diag(H))
    x <- x0
    for (j in (1 - S):(n * thinning)) {
      for (i in 1:d) {
        mu_i <- mean[i] - (1 / H[i, i]) * H[i,-i] %*% (x[-i] -
                                                         mean[-i])
        F_tmp <- pnorm(c(lower[i], upper[i]), mu_i, sd[i])
        Fa <- F_tmp[1]
        Fb <- F_tmp[2]
        x[i] <- mu_i + sd[i] * qnorm(U[l] * (Fb - Fa) + Fa)
        l <- l + 1
      }
      if (j > 0) {
        if (thinning == 1) {
          X[j,] <- x
        }
        else if (j %% thinning == 0) {
          X[j %/% thinning,] <- x
        }
      }
    }
    return(X)
  }



#' Simulate Probit Gaussian Process
#' @description Simulate Gaussian process replicates. A corresponding binary
#' random field is constructed by constructing an indicator for whether the
#' Gaussian process is > 0 at each location.
#' @param nrep Number of replicate processes to simulate
#' @param gvar Matern covariance function variance parameter
#' @param gscl Matern covariance function scale parameter
#' @param nu Matern smoothness parameter
#' @param nloc (optional) number of locations to simulate. Must be passed if
#' \code{locs} matrix is not passed.
#' @param locs (optional) nloc x 2 matrix of spatial locations at which to
#' simulate process
#' @param mu mean of Gaussian process. Should be either a (1) scalar, (2) vector
#' of length nloc, or (3) matrix of dimension nloc x nrep
#'
#' @return \code{nrep} replicates of Gaussian Process and 0 thresholded
#' Gaussian process matrix. Rows: different locations. Columns: independent
#' replicates
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(cowplot)
#' gp_list <- rbingp(nrep = 2, gvar = 1, gscl = 0.1, mu = -1, nloc = 1000)
#' bin_gp <- gp_list$bin_gp
#' bin <- gp_list$bin
#' s <- gp_list$locs
#' # Plot --------------------------------------------------------------------
#' p1 <- qplot(s1, s2, data=data.frame(s1 = s[,1], s2 = s[,2]), colour=bin_gp[,1]) +
#'   scale_color_gradient(low = "blue", high = "red")
#' p2 <- qplot(s1, s2, data=data.frame(s1 = s[,1], s2 = s[,2]), colour=bin[,1]) +
#'   scale_color_gradient(low = "blue", high = "red")
#' cowplot::plot_grid(p1, p2)
dgp   <-  function(Z,
                   locs,
                   gvar,
                   gscl,
                   nu = 1 / 2,
                   mu = 0) {
  RandomFields::RFoptions(spConform = FALSE)
  model <-
    RandomFields::RMmatern(
      var = gvar,
      scale = gscl,
      nu = nu,
      notinvnu = TRUE
    )
  return(RandomFields::RFlikelihood(model, x = locs, data = Z)$loglikelihood)
}

#' Calculate multivariate normal density frhom cholesky decomp. of Sigma
#' @description This function is an abbreviated form of the dmvnorm function
#' of the mvtnorm package when the cholesky factor has already been computed.
#' @param x matrix of multivariate normal observations. Indep. obs in each row
#' @param mean (vector) multivariate normal mean.
#' @param U (upper triangular) Cholesky decomposition of covariance matrix.
#' @param log (logical) should log be returned?
#' @export
#' @examples
#' library(mvtnorm)
#' n <- 400
#' A <- matrix(runif(n*n), n,n)
#' S <- (t(A) + A)/2 + diag(n)*n
#' H <- round(solve(S), 10)
#' mu <- x
#' x <- rmvnorm(5, mean = mu, sigma = S)
#' U <- chol(S)
#' ll1 <- dmvnorm(x, mu, S, log = TRUE)
#' ll2 <- dmvnormc(x, mu, U, log= TRUE)
#' all.equal(ll1,ll2)
dmvnorm_chol <-
  function (x,
            mean = rep(0, p),
            U = diag(p),
            log = FALSE)
  {
    if (is.vector(x))
      x <- matrix(x, ncol = length(x))
    p <- ncol(x)
    if (!missing(mean)) {
      if (!is.null(dim(mean)))
        dim(mean) <- NULL
      if (length(mean) != p)
        stop("mean and sigma have non-conforming size")
    }
    if (!missing(U)) {
      if (p != ncol(U))
        stop("x and sigma have non-conforming size")
      U[lower.tri(U)] <- 0    # ensure upper triangular
    }
    if (inherits(U, "error")) {
      x.is.mu <- colSums(t(x) != mean) == 0
      logretval <- rep.int(-Inf, nrow(x))
      logretval[x.is.mu] <- Inf
    }
    else {
      tmp <- backsolve(U, t(x) - mean, transpose = TRUE)
      rss <- colSums(tmp ^ 2)
      logretval <- -sum(log(diag(U))) - 0.5 * p * log(2 *
                                                        pi) - 0.5 * rss
    }
    names(logretval) <- rownames(x)
    if (log)
      logretval
    else
      exp(logretval)
  }


#' Inverse Gamma density
#' @description Note that the scale parameter of the inverse gamma is the rate
#' parameter of the corresponding gamma distribution. If W = 1/X, for
#' X ~ Gamma(shape= a, scale = b), then W ~ Inv. Gamma(shape = a, scale = 1/b)
#' @param x (scalar or vector) argument to inverse gamma density
#' @param shape parameter
#' @param rate 1/scale
#' @param scale 1/rate
#' @param log (logical) should log-density be returned
#' @return inverse gamma density at x.
#' @export
#'
#' @examples
#' a <- 9
#' b <- 7
#' curve(dinvgamma(x, shape = a, scale = b), from = 0, to = 5, col = 2)
#' y<- 1/rgamma(100000, shape = a, rate = b)
#' lines(density(y))
dinvgamma <-
  function (x,
            shape,
            rate = 1,
            scale = 1 / rate,
            log = FALSE)
  {
    if (missing(rate) && !missing(scale))
      rate <- 1 / scale
    log_f <- dgamma(1 / x, shape, scale = rate, log = TRUE) - 2 * log(x)
    if (log)
      return(log_f)
    exp(log_f)
  }

#' Calculate multivariate normal density from precision matrix
#'
#' @param x (vector or matrix) argument of MVN density. If a matrix, independent
#'          replicates belong to separate rows.
#' @param mu (vector) MVN mean
#' @param Omega (positive definite matrix) MVN precision
#' @param logdetOmega (optional scalar) log-determinant of precision matrix
#' @param log (logical) should log-likelihood be returned? Default = FALSE
#'
#' @return multivariate normal densities
#' @export
#'
#' @examples
#' library(mvtnorm)
#' n <- 100
#' p <- 2
#' A <- matrix(runif(p^2), nrow = p, ncol= p)
#' A <- (t(A) + A)/2 + diag(p)*p
#' iA <- solve(A)
#' mu <- rnorm(p)
#' x <- rmvnorm(n, mean = mu, sigma = A)
#' expected = dmvnorm(x, mu, A, log = T)
#' value = dmvnp(x, mu, iA, log = T)
dmvnp <- function (x, mu, Omega, logdetOmega = NULL, log = FALSE)
{
  if (!is.matrix(x))
    x <- rbind(x)
  if (!is.matrix(mu))
    mu <- matrix(mu, nrow(x), ncol(x), byrow = TRUE)
  if (missing(Omega))
    Omega <- diag(ncol(x))
  if (!is.matrix(Omega))
    Omega <- matrix(Omega)
  # if (!is.positive.definite(Omega))
  #   stop("Matrix Omega is not positive-definite.")
  k <- nrow(Omega)
  if(is.null(logdetOmega)){
    logdetOmega <- logdet(Omega)
  }
  ss <- x - mu
  z <- rowSums({
    ss %*% Omega
  } * ss)
  dens <- as.vector((-k/2) * log(2 * pi) + 0.5 * logdetOmega -
                      0.5 * z)
  if (log == FALSE)
    dens <- exp(dens)
  return(dens)
}

#' @export
logdet <- function(x){
  return(2 * sum(log(diag(cholC(x)))))
}



#' array index to single value index of matrix
#'
#' @param r row number
#' @param c column number
#' @param m total number of rows in matrix
#'
#' @return single number index corresponding to c,r
#' @export
sub2ind <- function(r,c, m){
  ind = (c-1)*m + r
}

#' Vectorized array index to single value conversion
#'
#' @param rs vector of row indices
#' @param cs vector of column indices
#' @param m total number of columns in matrix
#' @return vector of same length as rs and cs giving flattened matrix indices
#' that match those of the row, column pairs for c(rs[1], cs[1]),
#' c(rs[2], cs[2]), ... etc.
#' @export
#'
#' @examples
#' A <- matrix(1:9, nrow = 3, ncol = 3)
#' A[vsub2ind(1:3, c(3,1,2), 3)]     # should give 7,2, 6
vsub2ind <- function(rs, cs, m) {
  mapply(sub2ind,
         r = rs,
         c = cs,
         MoreArgs = list(m = m))
}

#' Transform from the range parameter to the microergodic parameter
#' @param range_par range parameter of matern
#'
#' @param smoothness_par smoothness parameter of matern
#' @param var_par variance parameter of matern
#'
#' @export
range_to_microerg <-
  function(range_par, smoothness_par, var_par = 1) {
    return(var_par / (range_par ^ (2 * smoothness_par)))
  }


#' Transform from the microergodic parameter to the smoothness parameter in matern
#' @param microerg_par microergodic parameter (See Kauffman and Shaby, 2012)
#'
#' @param smoothness_par smoothness parameter of matern
#' @param var_par variance parameter of matern
#'
#' @export
microerg_to_smooth <-
  function(microerg_par, range_par, var_par = 1) {
    return(log(var_par/microerg_par)/(2*log(range_par)))
  }


#' Calculate pairwise distances between rows of A matrix and rows of B matrix
#' @export calc_dist
#' @param A Matrix 1 (nrow = m)
#' @param B Matrix 2 (nrow = n)
#' @return An m by n matrix of distances, where m = \code{nrow(x)}
#'   and n = \code{nrow(y)}
# calc_dist <- function(x, y) {
#   x <- as.matrix(x)
#   y <- as.matrix(y)
#   m <- nrow(x)
#   n <- nrow(y)
#   xn = apply(x, 1, function(rvec)
#     crossprod(rvec, rvec))
#   yn = apply(y, 1, function(rvec)
#     crossprod(rvec, rvec))
#   tmp = matrix(rep(xn, n), nrow = m)
#   tmp = tmp +  matrix(rep(yn, m), nrow = m, byrow = TRUE)
#   sqrt(tmp - 2 * tcrossprod(x, y))
# }
calc_distR <- function(A,B) {
  # A: matrix with obersvation vectors
  #         (nrow = number of observations)
  #
  # B: matrix with another set of vectors
  #          (e.g. cluster centers)
  A <- as.matrix(A)
  B <- as.matrix(B)
  result = matrix(ncol=nrow(B), nrow=nrow(A))
  for (i in 1:nrow(A))
    for (j in 1:nrow(B))
      result[i,j] = sqrt(sum( (A[i,] - B[j,])^2 ))

  result
}

#' Thinplate spline kernel function
#'
#' @param d distance between x and knot location
#'
#' @return thin-plate kernel function evaluated for distance d
#' @export
thinplate_kernel <- function(d) {
  return((d ^ 2) * log(d))
}

#' Make basis matrix
#'
#' @param obs_coord n x 2 matrix of observation coordinates
#' @param knot_coord p x 2 matrix of knot locations
#' @param type string describing type of kernel to use: one of "thinplate" or
#' "gauss"
#' @param kern_bw (float) kernel bandwidth to use
#'
#' @export
#'
#' @examples
#' knot_coord <- seq(0, 1, l = 10)
#' knot_coord <- as.matrix(expand.grid(knot_coord, knot_coord)) + 1e-5
#' obs_coord <- seq(0, 1, l = 50)
#' obs_coord <- as.matrix(expand.grid(obs_coord, obs_coord))
#' B <- make_basis(obs_coord, knot_coord, "thinplate")
#' y <- B%*%rnorm(nrow(knot_coord))
#' p <- raster_occur_gp(y, coord = obs_coord, col_lim = NULL)
#' p
make_basis <-
  function(obs_coord,
           knot_coord,
           type = c("thinplate", "gauss"),
           kern_bw = NULL) {
    type <- match.arg(type)
    dists <- calc_dist(obs_coord, knot_coord)
    if (type == "gauss") {
      B <- dnorm(dists, mean = 0, sd = kern_bw)
    }
    if (type == "thinplate") {
      B <- thinplate_kernel(dists)
    }
    return(B)
  }


#' @export
calc_procrustes <-
  function(X,
           which_times,
           traj_len,
           par_opt = NULL,
           n_core = NULL) {
    ntimes <- length(which_times)
    # D <- matrix(NA_real_, ntimes, ntimes)
    # Columns of and A and B should be different time lags
    # Rows of ... should be different loading coefficient dimensions
    if (is.null(n_core)) {
      n_core <- parallel::detectCores()
    }
    doMC::registerDoMC(cores = n_core)
    if(is.null(par_opt)){
      chunk_size <- ceiling(ntimes / n_core)
      par_opt <- list(chunkSize = chunk_size)
    }
    D <- foreach::foreach(i = which_times,
                          .combine = 'rbind',
                          .options.mpi = par_opt) %dopar% {
                            d <- rep(NA_real_, ntimes)
                            k <- 1
                            for (j in which_times) {
                              d[k] <-
                                pracma::procrustes(X[(i - traj_len + 1):i,],
                                                   X[(j - traj_len + 1):j,])$d
                              k <- k + 1
                            }
                            return(d)
                          }
    return(D)
  }



#' Add elements of sup_list to sub_list
#' @param sup_list high-priority list (will overwrite sub_list values)
#' @param sub_list low-priority list (will be overwritten by sup_list values)
#'
#' @return combine sub_list and sup_list. Any element that is in both sup_list
#' and sub_list will be overwritten by sup_list
#' @export
#'
#' @examples
#' const <- list(x = 1)
#' my_const <- list(y = 3)
#' update_list(my_const, const)
update_list <- function(sup_list, sub_list){
  add_names <- names(sup_list)
  if(!is.null(add_names)){
    for(name in add_names){
      sub_list[[name]] <- sup_list[[name]]
    }
  }
  return(sub_list)
}

#' @export
call_to_list <- function(list_or_call){
  if(class(list_or_call) == "call"){
    return(eval(list_or_call))
  }
}


#' Taken embedding of the rows of X.
#' @description The first row of the output will be nstep*ntime + 1
#' @param X Matrix to construct taken embedding from. The rows of \code{X} are
#' vectors of variables observed for each time point. The taken embedding stacks
#' lagged rows together horizontally to construct trajectories.
#'
#' @param nstep Number of time steps to take in constructing a trajectory.
#' @param stride size of steps backwards in time to take when constructing
#' lagged time trajectories of the rows of X.
#' @return Lagged matrix X. The original covariates for the kth row will be
#' in the rightmost columns.
#' @export
embed_taken <- function(X, nstep, stride) {
  start_pos <- nstep * stride + 1
  end_pos <- nrow(X)
  ret_X <- X[start_pos:end_pos, , drop = FALSE]
  for (step in 1:nstep) {
    ret_X <-
      cbind(X[(start_pos - step * stride):(end_pos - step * stride), ], ret_X)
  }
  return(ret_X)
}

#' Make analogue similarity matrix
#' @param dist_mat distance matrix. e.g. the i,j th element gives the distance
#'                 between the ith day and the jth day's atmos vars
#' @param gauss_var Variance to be used in un-normalized gaussian kernel
#' @param threshold (logical) Should the similarity matrix be thresholded to only include similarity
#'                  of very near neighbors (reduce bias of unrelated neighbors)?
#' @param top_n_neigh If threshold == TRUE, this specifies the number of nearest
#'                  neighbors to retain.
#' @return Similarity matrix where the columns are sum normalized to 1 and elements in
#'         a column j give similarity of others to that column.
#'
#' @export
make_analogue_matrix <- function(dist_mat,
                                 gauss_var = var(c(dist_mat)),
                                 threshold = TRUE,
                                 top_n_neigh = 5) {
  K <- gauss_kernel(dist_mat, theta = gauss_var)
  diag(K) <- 0
  if (threshold) {
    for (i in 1:ncol(dist_mat)) {
      set_to_zero <- order(K[, i], decreasing = T)[-c(1:top_n_neigh)]
      K[set_to_zero, i] <- 0
    }
  }
  col_sums <- apply(K, 2, sum)
  analog_sim_mat <- t(t(K) / col_sums)
  return(analog_sim_mat)
}



