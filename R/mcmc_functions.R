#' @title Automatic proposal variance tuning update
#' @export
update_var <- function(cur_var, acpt_rt, opt_rt, gamma1) {
  exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt))
}


#' Simulate from full conditional for beta recursive prior
#' @param B_cur matrix of beta vectors (rows: covariates (p), columns: time)
#' @param D similarity matrix (dimension n-time x n-time) the t_th column
#'          describes the similarity between the t-time and all other times.
#'          t_th column of D (excluding t_th element). Diagonal should be 0.
#'          Columns should sum to 1.
#'           is used to calculate weights for B_cur[,-t] matrixd columns.
#' @param Z Gaussian process observation matrix (rows: locations, columns:
#'          independent replicates (time))
#' @param X Matrix or array of covariates for mean function of Z
#'          (i.e. mu_Z = X %*%beta). If X is a matrix, the covariates are
#'          common for all replicates and the dimensions are nloc x ncovariates.
#'          If X is an array, the dimensions are nloc x ncovariates x nreplicate
#'          (i.e. the last dimension corresponds to the time replicate).
#' @param Z_cov nloc x nloc spatial covariance matrix for Gaussian process
#' @param prior_var prior variance for beta vector (assumed common for all
#'                  dimensions of beta)
#' @param ndraw (default = 1) number of full conditional draws to make
#'
#' @return draws from full conditional for beta. If multiple draws are made,
#' independent replicates are stored in each row (nrep x p)
#' @export
r_fc_beta_occur <- function(B_cur, D, Z, X, Z_cov, prior_var, use_recursive_prior = TRUE) {
  B_prop <- B_cur
  p <- dim(X)[2]
  nrep <- dim(Z)[2]
  Siginv <- solveC(Z_cov)
  XsigX <- matrix(0, p, p)
  XsigZ <- matrix(0, p, 1)
  if (!length(dim(X)) == 3) {
    stop("X should be a 3-dimensional array")
  }
  if(use_recursive_prior){
    for (i in 1:nrep) {
      XsigX <- t(X[, , i]) %*% Siginv %*% X[, , i]
      XsigZ <- t(X[, , i]) %*% Siginv %*% Z[, i]
      Bd <- (B_cur[, -i] %*% D[-i, i]) / prior_var
      beta_fc_cov <- solveC((XsigX + diag(p) * (1 / prior_var)))
      beta_fc_mean <- beta_fc_cov %*% (XsigZ + Bd)
      B_prop[, i] <-
        t(mvtnorm::rmvnorm(1, mean = beta_fc_mean, sigma = beta_fc_cov))
    }
  }
  else{
    for (i in 1:nrep) {
      XsigX <- t(X[, , i]) %*% Siginv %*% X[, , i]
      XsigZ <- t(X[, , i]) %*% Siginv %*% Z[, i]
      beta_fc_cov <- solveC((XsigX + diag(p) * (1 / prior_var)))
      beta_fc_mean <- beta_fc_cov %*% XsigZ
      B_prop[, i] <-
        t(mvtnorm::rmvnorm(1, mean = beta_fc_mean, sigma = beta_fc_cov))
    }
  }
  return(B_prop)
}

#' @export
rtmvtnorm <-
  function(mu,
           H,
           lower,
           upper,
           burn_in_samples = 0,
           Z_cur = NULL,
           initialize_x0 = FALSE,
           parallelize = TRUE,
           n_core = NULL,
           par_opt = NULL) {
    nrep <- ncol(mu)
    nloc <- nrow(mu)
    if (parallelize) {
      if (is.null(n_core)) {
        stop("Must pass number of cores to use if parallelizing")
      }
      doMC::registerDoMC(cores = n_core)
      if (initialize_x0) {
        Z_cur <- foreach::foreach(j = 1:nrep,
                                  .combine = 'cbind',
                                  .options.mpi = par_opt)     %dopar% {
                                    t(
                                      rtmvtnormC(
                                        1,
                                        mean = mu[, j],
                                        H = H,
                                        lower = lower[, j],
                                        upper = upper[, j],
                                        init = c(0),
                                        burn_in_samples = burn_in_samples,
                                        thinning = 1,
                                        initialize_x0 = T
                                      )
                                    )
                                  }
      }
      else{
        Z_cur <- foreach::foreach(j = 1:nrep,
                                  .combine = 'cbind',
                                  .options.mpi = par_opt) %dopar% {
                                    t(
                                      rtmvtnormC(
                                        1,
                                        mean = mu[, j],
                                        H = H,
                                        lower = lower[, j],
                                        upper = upper[, j],
                                        init = Z_cur[, j],
                                        burn_in_samples = burn_in_samples,
                                        thinning = 1,
                                        initialize_x0 = F
                                      )
                                    )
                                  }
      }
    }
    else{
      if (initialize_x0) {
        Z_cur <- matrix(NA_real_, nrow = nloc, ncol = nrep)
        for (j in 1:nrep) {
          Z_cur[, j] <-  t(
            rtmvtnormC(
              1,
              mean = mu[, j],
              H = H,
              lower = lower[, j],
              upper = upper[, j],
              init = c(0),
              burn_in_samples = burn_in_samples,
              thinning = 1,
              initialize_x0 = T
            )
          )
        }
      }
      else{
        for (j in 1:nrep) {
          Z_cur[, j] <- t(
            rtmvtnormC(
              1,
              mean = mu[, j],
              H = H,
              lower = lower[, j],
              upper = upper[, j],
              init = Z_cur[, j],
              burn_in_samples = burn_in_samples,
              thinning = 1,
              initialize_x0 = F
            )
          )
        }
      }
    }
    return(Z_cur)
  }


# Update scale parameter
#' @export
update_occur_scale_smooth <-
  function(scale_cur,
           prop_var,
           cov_model_cur,
           coord,
           Z_cur,
           smooth_cur,
           gp_scale_smooth_U) {
    prop <-
      exp(chol_prop(rbind(log(scale_cur), log(smooth_cur)), U = gp_scale_smooth_U, prop_var = prop_var))
    if (all(prop > 0)) {
      scale_prop <- prop[1]
      smooth_prop <- prop[2]
      cov_model_prop <-
        RandomFields::RMmatern(
          var = 1,
          scale = scale_prop,
          nu = smooth_prop,
          notinvnu = TRUE
        )
      ll_cur <- RandomFields::RFlikelihood(cov_model_cur,
                                           x = coord, data = Z_cur)$loglikelihood
      ll_prop <- RandomFields::RFlikelihood(cov_model_prop,
                                            x = coord, data = Z_cur)$loglikelihood
      if (log(runif(1)) < ll_prop - ll_cur) {
        scale_cur <- scale_prop
        smooth_cur <- smooth_prop
        cov_model_cur <- cov_model_prop
        acpt <- 1
      }
      else{
        acpt <- 0
      }
    }
    else{
      acpt <- 0
    }
    return(list(
      scale = scale_cur,
      smooth = smooth_cur,
      cov_model_cur = cov_model_cur,
      acpt = acpt
    ))
  }

#' Simulate from full conditional for beta (recursive prior)
#' @param B_cur matrix of beta vectors (rows: covariates (p), columns: time)
#' @param D similarity matrix (dimension n-time x n-time) the t_th column
#'          describes the similarity between the t-time and all other times.
#'          t_th column of D (excluding t_th element). Diagonal should be 0.
#'          Columns should sum to 1.
#'           is used to calculate weights for B_cur[,-t] matrixd columns.
#' @param Y_list list of (nloc x nrep) matrices of observations. List is
#'        of length nrep.
#' @param X_list List of covariates for mean function of Y. Each element of the
#'          list is for a different year. (i.e. mu_Y_t = X_t(s) %*%beta_t).
#'          The dimensions are nloc x p
#' @param gp_inv_cor_list List of (nloc x nloc) GP inverse correlation
#'        correlation matrices. List is of length nrep.
#' @param sigma_sq (vector of length = nrep) containing the random GP variance
#'        for each replicate.
#' @param prior_var prior variance for beta vector (assumed common for all
#'                  dimensions of beta)
#' @return draws from full conditional for beta. If multiple draws are made,
#' independent replicates are stored in each row (nrep x p)
#' @export
r_fc_beta_intensity <-
  function(B_cur,
           D,
           Y_list,
           X_list,
           clust_loc,
           gp_inv_cor_list,
           sigma_sq,
           prior_var,
           use_recursive_prior = TRUE) {
    B_prop <- B_cur
    p <- nrow(B_cur)
    nrep <- ncol(B_cur)
    if(use_recursive_prior){
      for (i in 1:nrep) {
        Y_list[[i]] <-  Y_list[[i]] - clust_loc[i]
        Siginv <- gp_inv_cor_list[[i]] / sigma_sq[i]
        XsigX <- t(X_list[[i]]) %*% Siginv %*% X_list[[i]]
        XsigY <- t(X_list[[i]]) %*% Siginv %*% Y_list[[i]]
        Bd <- (B_cur[, -i] %*% D[-i, i]) / prior_var
        beta_fc_cov <- solveC((XsigX + diag(p) * (1 / prior_var)))
        beta_fc_mean <- beta_fc_cov %*% (XsigY + Bd)
        B_prop[, i] <-
          t(mvtnorm::rmvnorm(1, mean = beta_fc_mean, sigma = beta_fc_cov))
      }
    }
    else{
      for (i in 1:nrep) {
        Y_list[[i]] <-  Y_list[[i]] - clust_loc[i]
        Siginv <- gp_inv_cor_list[[i]] / sigma_sq[i]
        XsigX <- t(X_list[[i]]) %*% Siginv %*% X_list[[i]]
        XsigY <- t(X_list[[i]]) %*% Siginv %*% Y_list[[i]]
        beta_fc_cov <- solveC((XsigX + diag(p) * (1 / prior_var)))
        beta_fc_mean <- beta_fc_cov %*% XsigY
        B_prop[, i] <-
          t(mvtnorm::rmvnorm(1, mean = beta_fc_mean, sigma = beta_fc_cov))
      }
    }
    return(B_prop)
  }

#' @export
r_fc_beta_prior_var <-
  function(beta_prior_var,
           B_cur,
           D,
           prior_ig_shape = 0.001,
           prior_ig_scale = 0.001,      # some shrinkage
           use_recursive_prior = TRUE) {
    p <- nrow(B_cur)
    nrep <- ncol(B_cur)
    if (use_recursive_prior) {
      ss <- 0
      for (i in 1:nrep) {
        ss <- ss + sum((B_cur[, i] - B_cur[,-i] %*% D[-i, i])^2)
      }
    }
    else{ # Indep. prior
      ss <- sum(B_cur^2)
    }
    return(1 / rgamma(1,
                      shape = prior_ig_shape + p *nrep / 2,
                      rate =  prior_ig_scale + ss / 2))
  }


#' @export
r_fc_sigma_sq_intensity <-
  function(Y_list, clust_loc, X_list, B, gp_inv_cor_list, a, b) {
    nrep <- ncol(B)
    nlocs <- sapply(Y_list, length)
    sigma_sq <- rep(NA_real_, nrep)
    alpha <- 1:nrep
    beta <- 1:nrep
    for (i in 1:nrep) {
      r <- Y_list[[i]] - clust_loc[i] - X_list[[i]] %*% B[, i]
      alpha[i] <- (a[i] + nlocs[i]) / 2
      beta[i] <- (a[i] * b[i] + t(r) %*% gp_inv_cor_list[[i]] %*% r) / 2
      sigma_sq[i] <- 1 / rgamma(1, shape = alpha[i], rate = beta[i])
    }
    return(sigma_sq)
  }

#' Adaptive metropolis update for scale
#' @export
update_scale_smooth_intensity <-
  function(scale_tilde_cur,
           prop_var,
           Y_list,
           clust_loc,
           X_list,
           B,
           sigma_sq,
           smooth_tilde_cur,
           gp_inv_cor_list_cur,
           coord,
           subset_id_list,
           clust_labels,
           unq_clust,
           scale_range,
           smooth_range,
           gp_scale_smooth_U) {
    nclust <- length(unq_clust)
    acpt <- rep(NA_real_, nclust)
    nrep <- ncol(B)
    smooth_tilde_prop <- scale_tilde_prop <- rep(NA_real_, nclust)
    gp_inv_cor_list_prop <- gp_inv_cor_list_cur
    for (k in 1:nclust) {
      this_clust <- which(clust_labels == k)
      prop <-
        exp(chol_prop(rbind(log(scale_tilde_cur[k]), log(smooth_tilde_cur[k])),
                      U = gp_scale_smooth_U[[k]], prop_var = prop_var[k]))
      scale_tilde_prop[k] <- prop[1]
      smooth_tilde_prop[k] <- prop[2]
      if ((scale_tilde_prop[k] > scale_range[1])&
          (scale_tilde_prop[k] < scale_range[2])&
          (smooth_tilde_prop[k] > smooth_range[1])&
          (smooth_tilde_prop[k] < smooth_range[2])) {
        if (any(this_clust)) {
          ll_prop <- 0
          ll_cur <- 0
          prop_model <-
            RandomFields::RMmatern(
              var = 1,
              scale = scale_tilde_prop[k],
              nu = smooth_tilde_prop[k],
              notinvnu = TRUE
            )
          gp_cor_tilde_prop <- RandomFields::RFcovmatrix(prop_model, x = coord)
          for (i in this_clust) {
            ll_prop <-
              tryCatch({
                gp_inv_cor_list_prop[[i]] <-
                  solveC(gp_cor_tilde_prop[subset_id_list[[i]],
                                           subset_id_list[[i]]])
                ll_prop <-
                  ll_prop + dmvnp(
                    c(
                      Y_list[[i]] - X_list[[i]] %*% B[, i] - clust_loc[i]
                    ),
                    mu =  rep(0, length(Y_list[[i]])),
                    Omega = gp_inv_cor_list_prop[[i]] / sigma_sq[i],
                    log = TRUE
                  )
              },
              error = function(e) {
                -Inf
              })

            ll_cur <-
              ll_cur + dmvnp(
                c(
                  Y_list[[i]] - X_list[[i]] %*% B[, i] - clust_loc[i]
                ),
                mu =  rep(0, length(Y_list[[i]])),
                Omega = gp_inv_cor_list_cur[[i]] /
                  sigma_sq[i],
                log = TRUE
              )
          }
          if (log(runif(1)) < ll_prop - ll_cur) {
            scale_tilde_cur[k] <- scale_tilde_prop[k]
            smooth_tilde_cur[k] <- smooth_tilde_prop[k]
            for (i in this_clust) {
              gp_inv_cor_list_cur[[i]] <- gp_inv_cor_list_prop[[i]]
            }
            acpt[k] <- 1
          }
          else{
            acpt[k] <- 0
          }
        }
        else{
          scale_tilde_cur[k] <- runif(1, scale_range[1], scale_range[2])
          smooth_tilde_cur[k] <- runif(1, smooth_range[1], smooth_range[2])
          acpt[k] <-  1                 # How to update acceptance vectors when cluster is empty?
        }
      }
      else{
        acpt[k] <- 0
      }
    }
    return(
      list(
        scale_tilde = scale_tilde_cur,
        smooth_tilde = smooth_tilde_cur,
        gp_inv_cor_list = gp_inv_cor_list_cur,
        acpt = acpt
      )
    )
  }

#' @export
lfc_a_b <-
  function(a,
           b,
           sigma_sq,
           prior_shape = 1,
           prior_scale = 1) {
    return(sum(dinvgamma(
      sigma_sq, a / 2, scale = a * b / 2, log = TRUE
    )) +
      dgamma(
        b,
        shape = prior_shape,
        scale = prior_scale,
        log = TRUE
      ))
  }

#' Update T d.f. parameter in mcmc
#' @export
#' @examples
#' nmcmc <- 1000
#' n <- 300
#' nclust <- 4
#' clust_labels <- sample(1:(nclust -1), size = n, replace = T)
#' a_tilde <- c(0, 2, 4, 5)/10 + 1
#' b_tilde <- rep(10, nclust)
#' a <- a_tilde[clust_labels]
#' b <- b_tilde[clust_labels]
#' sigma_sq <- 1 / rgamma(n, shape = a/2, rate = a*b/2)
#' a_up <- matrix(1, nmcmc, nclust)
#' for(i in 1:nmcmc){
#'   a_up[i,] <- update_a(a_tilde, b_tilde, sigma_sq, unq_clust = 1:nclust, clust_labels = clust_labels)
#' }
#' matplot(a_up, type = "l")
update_a <-
  function(a_tilde,
           b_tilde,
           sigma_sq,
           unq_clust,
           clust_labels,
           lower_bound = 0.2,
           upper_bound = 20,
           b_prior_shape = 1,
           b_prior_scale = 1) {
    nclust <- length(unq_clust)
    for (k in 1:nclust) {
      this_clust <- which(clust_labels == k)
      if (any(this_clust)) {
        a_prop <- rnorm(1, a_tilde[k], 0.1)
        if ((a_prop < lower_bound) | (a_prop > upper_bound)) {
          ll_prop <- -Inf
        }
        else{
          ll_prop <- lfc_a_b(a_prop, b_tilde[k], sigma_sq[this_clust],
                             prior_shape = b_prior_shape, prior_scale = b_prior_scale)
        }
        ll_cur <-
          lfc_a_b(a_tilde[k], b_tilde[k], sigma_sq[this_clust],
                  prior_shape = b_prior_shape, prior_scale = b_prior_scale)
        if (log(runif(1)) < ll_prop - ll_cur) {
          a_tilde[k] <- a_prop
        }
      }
      else{
        a_tilde[k] <-
          runif(1, lower_bound, upper_bound)   # uniform on (0,20) in 0.2 increments
      }
    }
    return(a_tilde)
  }

#' Function to update b_tilde in mcmc loop
#' @param b_tilde Vector of length nclust of scale parameters for
#' distributions for sigma_sq
#' @param a_tilde Vector of length nclust of scale parameters for
#' distributions for sigma_sq
#' @param sigma_sq Vector of length nrep of GP variances
#' @param unq_clust vector giving unique clusters (e.g. 1:nclust)
#' @param clust_labels vector of length nrep of cluster labels for each
#' replicate
#' @param prior_shape gamma hyperprior shape parameter
#' @param prior_scale gamma hyperprior scale parameter
#' @export
#' @examples
#' nmcmc <- 1000
#' n <- 300
#' nclust <- 4
#' clust_labels <- sample(1:nclust, size = n, replace = T)
#' a_tilde <- rep(0.1, nclust)
#' b_tilde <- c(1, 10, 20,35)
#' a <- a_tilde[clust_labels]
#' b <- b_tilde[clust_labels]
#' sigma_sq <- 1 / rgamma(n, shape = a/2, rate = a*b/2)
#' b_up <- matrix(1, nmcmc, nclust)
#' for(i in 1:nmcmc){
#'   b_up[i,] <- update_b(b_tilde,a_tilde, sigma_sq, 1:nclust,clust_labels,prior_shape = 1, prior_scale = 100)
#' }
#' matplot(b_up, type = "l")
#' apply(b_up,2, mean)
update_b <-
  function(b_tilde,
           a_tilde,
           sigma_sq,
           unq_clust,
           clust_labels,
           prior_shape = 1,
           prior_scale = 1) {
    nclust <- length(unq_clust)
    for (k in 1:nclust) {
      this_clust <- which(clust_labels == k)
      if (any(this_clust)) {
        alpha0 <- a_tilde[k] * length(this_clust) / 2 + prior_shape
        beta0 <-
          1 / ((a_tilde[k] / 2) * sum(1 / sigma_sq[this_clust]) + 1 / prior_scale)
        b_tilde[k] <- rgamma(1, alpha0, scale =  beta0)
      }
      else{
        b_tilde[k] <- rgamma(1, shape = prior_shape, scale = prior_scale)
      }
    }
    return(b_tilde)
  }

#' Multivariate normal random walk proposals on whitened space
#'
#' @param x p dimensional vector (current value in random walk)
#' @param U upper triangular matrix of Cholesky decomposition p x p
#' @param prop_var (scalar) may be adapted
#'
#' @return Multivariate normal random walk proposal of dimension p
#' @export
#' @examples
#' ### not run
#' library(mvtnorm)
#' set.seed(1)
#' n <- 2
#' S <- cbind(c(1,4/5), c(4/5, 1))
#' U <- chol(S)
#' z <- t(rmvnorm(100, sigma = S))
#' znew <- chol_prop(z, U, 0.0001)
#' plot(t(z))
#' points(t(znew), col = 2)
chol_prop <- function(x, U, prop_var) {
  return(t(U) %*% ((t(solveC(
    U
  )) %*% x) +
    rnorm(
      prod(dim(x)),
      mean = 0,
      sd = matrix(sqrt(prop_var), ncol = ncol(x), nrow = nrow(x), byrow = T)
    )))
}

#' Calculate class probabilites
#' @param X T x p matrix of covariates (T = number of replicates, p = number of
#' covariates).
#' @param A p x (K-1) matrix of covariate loadings. The Kth class is the
#' reference class

#' @export
calc_class_probs <- function(X, A) {
  numers <- cbind(exp(X %*% A), 1)
  return(numers / rowSums(numers))
}


#' Log-likelihood for categorical dist
#'
#' @param P_mat T x K matrix of probabilities (T = n-replicates, K = n-classes)
#' @param labels T vector of class labels each in 1:K
#'
#' @return categorical distribution log-likelihood
#' @export
#'
#' @examples
#' X <- matrix(c(1,0, 1, 0,1,1), nrow = 3, ncol = 2)  # p = 2, T = 3
#' A <- matrix(1:8, nrow = 2, ncol = 4) # p= 2, K = 5
#' P <- calc_class_probs(X ,A)
#' xi <- c(1,3, 5)             # class labels for replicates 1,2,3
#' l_categorical(P, xi)
l_categorical <- function(P_mat, labels) {
  R <- nrow(P_mat)
  return(sum(log(P_mat[vsub2ind(1:R, labels, R)])))
}

#' @export
l_fc_alpha <-
  function(A, X, clust_labels, alpha_prior_var, Pr = NULL) {
    if (is.null(Pr)) {
      Pr <- calc_class_probs(X, A)
    }
    return(l_categorical(Pr, clust_labels) +
             sum(dnorm(
               A, sd = sqrt(alpha_prior_var), log = T
             )))
  }

#' @export
update_alpha <-
  function(A_cur,
           X,
           clust_labels,
           prop_var_mat,
           alpha_prior_var) {
    nclust <- ncol(A_cur) + 1
    p <- nrow(A_cur)
    A_prop <- A_cur
    acpt <- matrix(0, nrow = p, ncol = nclust - 1)
    for (i in 1:length(A_cur)) {
      A_prop[i] <- rnorm(1, A_cur[i], sd = sqrt(prop_var_mat[i]))
      if (log(runif(1)) <
          l_fc_alpha(A_prop, X, clust_labels, alpha_prior_var) -
          l_fc_alpha(A_cur, X, clust_labels, alpha_prior_var)) {
        A_cur[i] <- A_prop[i]
        acpt[i] <- 1
      }
      else{
        A_prop[i] <- A_cur[i]
      }
    }
    return(list(A = A_cur,
                acpt = acpt))
  }


llik_one <-
  function(y,
           sigma_sq,
           clust_loc,
           mu,
           gp_inv_cor,
           a,
           b) {
    return(
        dmvnp(
          c(y - clust_loc - mu),
          mu =  rep(0, length(y)),
          Omega = gp_inv_cor / sigma_sq,
          log = TRUE
        ) +
        dinvgamma(
          sigma_sq,
          shape = a / 2,
          scale = a * b / 2,
          log = TRUE
        )
    )
  }

l_fc_clust_label_one <-
  function(clust_label,
           clust_probs,
           y,
           sigma_sq,
           clust_loc,
           mu,
           gp_inv_cor,
           a,
           b) {
    return(llik_one(y, sigma_sq, clust_loc, mu, gp_inv_cor, a, b) +
      log(clust_probs[clust_label]))
  }


#' @export
update_clust_labels <-
  function(clust_labels_cur,
           unq_clust,
           clust_probs,
           Y_list,
           clust_loc_tilde,
           sigma_sq,
           gp_inv_cor_list,
           a_tilde,
           b_tilde,
           X_list,
           B,
           gp_cor_tilde_list_cur,
           subset_id_list) {
    nrep <- length(Y_list)
    acpt <- rep(0, nrep)
    for(i in 1:nrep){
      clust_label_prop <- sample(unq_clust, 1)
      if(clust_label_prop != clust_labels_cur[i]){
        ll_prop <-
          tryCatch({
            gp_inv_cor_prop <-
              solveC(gp_cor_tilde_list_cur[[clust_label_prop]][subset_id_list[[i]],
                                                               subset_id_list[[i]]])
            ll_prop <- l_fc_clust_label_one(
              clust_label_prop,
              clust_probs[i, ],
              Y_list[[i]],
              sigma_sq[i],
              clust_loc_tilde[clust_label_prop],
              X_list[[i]] %*% B[, i],
              gp_inv_cor_prop,
              a_tilde[clust_label_prop],
              b_tilde[clust_label_prop]
            )
          },
          error = function(e) {
            -Inf
          })
      ll_cur <- l_fc_clust_label_one(
        clust_labels_cur[i],
        clust_probs[i, ],
        Y_list[[i]],
        sigma_sq[i],
        clust_loc_tilde[clust_labels_cur[i]],
        X_list[[i]] %*% B[, i],
        gp_inv_cor_list[[i]],
        a_tilde[clust_labels_cur[i]],
        b_tilde[clust_labels_cur[i]]
      )
      if(log(runif(1)) < ll_prop - ll_cur){
        clust_labels_cur[i] <- clust_label_prop
        acpt[i] <- 1
        gp_inv_cor_list[[i]] <- gp_inv_cor_prop
      }
      }
    }
    return(list(clust_labels_cur = clust_labels_cur,
                acpt = acpt,
                gp_inv_cor_list = gp_inv_cor_list))
  }



#' @export
r_fc_clust_loc_tilde <-
  function(clust_loc_tilde,
           Y_list,
           B,
           X_list,
           gp_inv_cor_list,
           sigma_sq,
           clust_labels,
           prior_var,
           prior_clust_loc_tilde_mean = rep(0, length(clust_loc_tilde))
           ) {
    clust_loc_tilde_prop <- clust_loc_tilde
    for (k in 1:length(clust_loc_tilde)) {
      this_clust <- which(clust_labels == k)
      iv <- 1 / prior_var
      num <- prior_clust_loc_tilde_mean[k]/prior_var
      if (any(this_clust)) {
        for (i in this_clust) {
          w <- Y_list[[i]] - X_list[[i]] %*% B[, i]
          num <-
            num + t(w) %*% (gp_inv_cor_list[[i]] / sigma_sq[i]) %*% rep(1, length(w))
          iv <- iv + sum((gp_inv_cor_list[[i]] / sigma_sq[i]))
        }
      }
      clust_loc_tilde_prop[k] <-
        rnorm(1, mean = num / iv, sd = sqrt(1 / iv))
    }
    return(clust_loc_tilde_prop)
  }

