#' TODO: Add updates for the smoothness parameter. Automatically detect if
#' the certain parameters are passed in \code{const} list. If they are, then
#' don't update them. Otherwise, update them and require that initial values
#' be set in \code{init}.
#' @param Y (nloc x nrep) matrix of 0/1s indicating absence/occurrence of an
#' event. This is the response in the Spatial GLM.
#' @param X (nloc x p x nrep) array of covariates where p is the number of
#' covariates.
#' @param D (nrep x nrep) matrix of similarity measures. The i,j element gives
#' the similarity between the jth replicate and the ith replicate. Each column
#' should sum to one.
#' @param coord (nloc x 2) matrix of observations coordinates
#' @param init list of initial values for the parameters.
#' Should include:
#'        beta, gp_scale
#' Optionally:
#'        Z matrix of latent GP (dimension matches Y)
#' @param tune_var list of proposal variance initial values.
#' Should include:
#'        gp_scale (GP range)
#' @param const list of constants.
#' @param niter (int) number of MCMC iterations
#' @param thin_int (int) Thinning interval
#' @param opt_rt (0 < opt_rt < 1) target acceptance rate for M-H updated
#' parameters. Used to tune proposal variances.
#' @param n_core Number of cores to parallelize over. If this is not set and
#' \code{parallelize = TRUE}, n_core will be set to the number of available
#' cores.
#' @param parallelize (logical) Should independent latent Z processes be
#' updated in parallel?
#' @param progress_bar (logical) Should a progress bar be printed to the
#' terminal? Good for short runs.
#'
#' @export
mcmc_occur <-
  function(Y,
           X,
           D,
           coord,
           init,
           tune_var = list(gp_scale_smooth = 0.1),
           const = list(beta_prior_var = 1),
           niter = 100,
           thin_int = 1,
           opt_rt = 0.4,
           n_core = NULL,
           parallelize = TRUE,
           progress_bar = FALSE,
           save_on_error = TRUE,
           use_recursive_prior = TRUE) {
    if(save_on_error){
      options(error = quote(dump.frames("occur_dump", to.file = TRUE, include.GlobalEnv = TRUE)))
    }

    # Update Defaults ---------------------------------------------------------
    # Get the default arguments in this function and set them
    form_arg <- formals(mcmc_occur)
    const <- update_list(const, form_arg$const)
    tune_var <- update_list(tune_var, form_arg$tune_var)

    # Setup parallelization
    if (parallelize) {
      if (is.null(n_core)) {
        n_core = parallel::detectCores()
      }
      chunk_size <- ceiling(nrep / n_core)
      par_opt <- list(chunkSize = chunk_size)
    }
    p <- dim(X)[2]
    nrep <- ncol(Y)
    nloc <- nrow(Y)
    nsamples <- floor(niter / thin_int)

    # Matrices to store mcmc samples
    samples <- list(
      beta = array(NA_real_, dim = c(nsamples, p, nrep)),
      gp_scale =  rep( NA_real_, nsamples),
      gp_smooth = rep(NA_real_, nsamples),
      Z = array(NA_real_, dim = c(nsamples, nloc, nrep)),
      beta_prior_var = rep(NA_real_, nsamples)
    )

    # Automatic Tuning of Proposal variance
    c0 <- 10
    c1 <- 0.8
    tune_k <- 3
    win_len <-
      min(niter, 50)                       # Window for acceptance rates and prop vars
    # Proposal Variances
    prop_var <- list(
      gp_scale_smooth = c(tune_var$gp_scale_smooth, rep(NA_real_, win_len - 1))
    )
    acpt <- list(
      gp_scale = rep(NA_real_, win_len)
    )
    # Initialize cholesky list
    gp_scale_smooth_U <- diag(2)

    # Get initial values
    mu_cur <- matrix(NA_real_, nloc, nrep)
    beta_cur <- init$beta
    gp_scale_cur <- init$gp_scale
    gp_smooth_cur <- init$gp_smooth
    # If an initial value is passed for the beta_prior_var, then sample, otherwise fix
    beta_prior_var_cur <-ifelse(is.null(init$beta_prior_var),
                                const$beta_prior_var, init$beta_prior_var)

    # Latent process
    Z_cur <- matrix(NA_real_, nrow = nloc, ncol = nrep)
    Z_lb <- ifelse(Y > 0, 0, -Inf)
    Z_ub <- ifelse(Y == 0, 0, Inf)
    cov_model_cur <-
      RandomFields::RMmatern(
        var = 1,
        scale = gp_scale_cur,
        nu = gp_smooth_cur,
        notinvnu = TRUE
      )
    Z_cov_cur <-
      RandomFields::RFcovmatrix(model = cov_model_cur, x = coord)
    H_cov_cur <- solve(Z_cov_cur)
    for (j in 1:nrep) {
      mu_cur[, j] <- X[, , j] %*% beta_cur[, j]
    }
    if (!is.null(init$Z)) {
      Z_cur <- init$Z
    }
    else{
      Z_cur <-
        rtmvtnorm(
          mu_cur,
          H_cov_cur,
          Z_lb,
          Z_ub,
          burn_in_samples = 100,
          initialize_x0 = TRUE,
          parallelize = parallelize,
          n_core = n_core,
          par_opt = par_opt
        )
    }
    ### Begin MCMC loop ###
    smp_loc <- 1
    if (progress_bar) {
      pb <- txtProgressBar(min = 0,
                           max = nsamples,
                           style = 3)
    }
    for (i in 1:niter) {
      gamma1 <-
        c0 / (i + tune_k) ^ (c1)                # Automatic tuning constant
      if (progress_bar) {
        setTxtProgressBar(pb, smp_loc)
      }
      pos <- (i - 1) %% win_len + 1
      nxt_pos <- i %% win_len + 1

      # Calculate current mean -------------------------------------------------
      for (j in 1:nrep) {
        mu_cur[, j] <- X[, , j] %*% beta_cur[, j]
      }
      # Update latent Gaussian process -----------------------------------------
      Z_cur <-
        rtmvtnorm(
          mu_cur,
          H_cov_cur,
          Z_lb,
          Z_ub,
          burn_in_samples = 2,
          Z_cur = Z_cur,
          initialize_x0 = FALSE,
          parallelize = parallelize,
          n_core = n_core,
          par_opt = par_opt
        )
      # Update Beta ------------------------------------------------------------
      beta_cur <-
        r_fc_beta_occur(beta_cur,
                  D,
                  Z_cur,
                  X,
                  Z_cov = Z_cov_cur,
                  prior_var = beta_prior_var_cur)

      # Update beta prior variance if init is passed ----------------------------
      if(!is.null(init$beta_prior_var)){
        beta_prior_var_cur <- r_fc_beta_prior_var(beta_prior_var = beta_prior_var_cur,
                                                  B_cur = beta_cur,
                                                  D = D,
                                                  use_recursive_prior = use_recursive_prior)
      }

      # Update scale and smooth parameters
      if (i > thin_int * 5) {
        gp_scale_smooth_cov <- cov(cbind(log(samples$gp_scale),
                                         log(samples$gp_smooth)),
                                   use = 'complete.obs') +  diag(2)
        gp_scale_smooth_U <- tryCatch(chol(gp_scale_smooth_cov),
                                      error = function(e)
                                        diag(2)
                                    )
      }

      gp_scale_smooth_updates <- update_occur_scale_smooth(
        scale_cur = gp_scale_cur,
        prop_var = prop_var$gp_scale_smooth[pos],
        cov_model_cur = cov_model_cur,
        coord = coord,
        Z_cur = Z_cur,
        smooth_cur = gp_smooth_cur,
        gp_scale_smooth_U = gp_scale_smooth_U
      )
      acpt$gp_scale[pos] <- gp_scale_smooth_updates$acpt
      prop_var$gp_scale_smooth[nxt_pos] <-
        update_var(
          prop_var$gp_scale_smooth[pos],
          acpt_rt = mean(acpt$gp_scale, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )
      if (gp_scale_cur != gp_scale_smooth_updates$scale) {
        gp_scale_cur <- gp_scale_smooth_updates$scale
        gp_smooth_cur <- gp_scale_smooth_updates$smooth
        cov_model_cur <- gp_scale_smooth_updates$cov_model_cur
        Z_cov_cur <-
          RandomFields::RFcovmatrix(cov_model_cur, x = coord)
        H_cov_cur <- solve(Z_cov_cur)
      }

      # Save samples------------------------------------------------------------
      if (i %% thin_int == 0) {
        samples$Z[smp_loc, ,] <- Z_cur
        samples$beta[smp_loc, , ] <- beta_cur
        samples$gp_scale[smp_loc] <- gp_scale_cur
        samples$gp_smooth[smp_loc] <- gp_smooth_cur
        samples$beta_prior_var[smp_loc] <-  beta_prior_var_cur
        smp_loc <- smp_loc  + 1
      }
    }
    acpt_rts = mean(acpt$gp_scale, na.rm = TRUE)
    return(
      list(
        samples = samples,
        prop_var = prop_var,
        acpt_rts = acpt_rts,
        const = const
      )
    )
  }
