#' @export
mcmc_intensity <-
  function(Y_list,
           X_list,
           X_class,      # matrix dim = T x p_class_covariates
           D,
           coord,
           subset_id_list,
           nclust,
           init,
           tune_var = list(gp_scale_smooth_tilde = 0.1,
                           sigma_sq = 0.1,
                           alpha = 0.1),
           const = list(alpha_prior_var = 10,
                        prior_clust_loc_tilde_mean = NULL,
                        clust_loc_prior_var = 1,
                        beta_prior_var = 1,
                        b_prior_scale = 0.1,
                        b_prior_shape = 0.1,
                        a_lower_bound = 10,
                        a_upper_bound = 20,
                        gp_scale_range = NULL,
                        gp_smooth_range = c(1/2,3/2)),
           niter = 100,
           thin_int = 1,
           opt_rt = 0.4,
           n_core = NULL,
           parallelize = TRUE,
           progress_bar = FALSE,
           save_on_error = TRUE,
           use_recursive_prior = TRUE) {
    if(save_on_error){
      options(error = quote(dump.frames("intensity_dump", to.file = TRUE, include.GlobalEnv = TRUE)))
      # to load, run load("intensity_dump.rda") and debugger(intensity_dump)
    }

    # Overwrite any defaults in lists that have been passed -------------------
    # Get the default arguments in this function and set them
    form_arg <- formals(mcmc_intensity)
    form_arg$const <- call_to_list(form_arg$const)
    form_arg$tune_var <- call_to_list(form_arg$tune_var)
    const <- update_list(const, form_arg$const)
    tune_var <- update_list(tune_var, form_arg$tune_var)
    if(is.null(const$gp_scale_range)){
      coord_dist <- dist(coord)
      const$gp_scale_range <- c(min(coord_dist[coord_dist>0])/2, 4*max(coord_dist))
    }
    if(is.null(const$prior_clust_loc_tilde_mean)){
      const$prior_clust_loc_tilde_mean <-
        unname(quantile(sapply(Y_list, mean, na.rm = T),
                        seq(1/nclust, 1-1/nclust, l = nclust)))
    }

    # Make sure that the similarity matrix is normalized ---------------------
    if(any(abs(colSums(D) - 1) > .Machine$double.eps^0.9)){
      warning("Columns of similarity matrix must be normalized to have sum = 1. Re-normalizing.")
      D <- D %*% diag((1/colSums(D)))
    }

    # Define constants---------------------------------------------------------
    p <- nrow(init$B)
    p_class_covariates <- ncol(X_class)
    nrep <- length(Y_list)
    nsamples <- floor(niter / thin_int)
    unq_clust = 1:nclust

    # Setup parallelization ----------------------------------------------------
    if (parallelize) {
      if (is.null(n_core)) {
        n_core = parallel::detectCores()
      }
      chunk_size <- ceiling(nrep / n_core)
      par_opt <- list(chunkSize = chunk_size)
    }

    # Matrices to store mcmc samples--------------------------------------------
    samples <- list(
      clust_labels = matrix(NA_real_, nrow = nsamples, ncol = nrep),
      a_tilde = matrix(NA_real_, nrow = nsamples, ncol = nclust),
      b_tilde = matrix(NA_real_, nrow = nsamples, ncol = nclust),
      gp_scale_tilde =  matrix( NA_real_, nrow = nsamples, ncol = nclust),
      gp_smooth_tilde = matrix(NA_real_, nrow = nsamples, ncol = nclust),
      clust_loc_tilde = matrix(NA_real_, nrow = nsamples, nclust),
      sigma_sq = matrix(NA_real_, nrow = nsamples, ncol = nrep),
      B = array(NA_real_, dim = c(nsamples, p, nrep)),
      alpha = array(NA_real_, dim = c(nsamples, p_class_covariates, nclust - 1)),
      beta_prior_var = rep(NA_real_, nsamples)
    )

    # Automatic Tuning of Proposal variance-------------------------------------
    c0 <- 10
    c1 <- 0.8
    tune_k <- 3
    win_len <- min(niter, 50)

    # Proposal Variances--------------------------------------------------------
    prop_var <- list(
      gp_scale_smooth_tilde = matrix(NA_real_, nrow = win_len, ncol = nclust),
      alpha = array(NA_real_, dim = c(win_len, p_class_covariates, nclust - 1))
    )
    prop_var$gp_scale_smooth_tilde[1,] <- tune_var$gp_scale_smooth_tilde
    prop_var$alpha[1,,] <- tune_var$alpha
    acpt <- list(
      gp_scale_smooth_tilde = matrix(NA_real_, nrow = win_len, ncol = nclust),
      sigma_sq = matrix(NA_real_, nrow = win_len, ncol = nrep),
      alpha = array(NA_real_, dim = c(win_len, p_class_covariates, nclust - 1)),
      clust_labels = matrix(NA_real_, win_len, nrep)
    )

    # Initialize cholesky list for scale/smooth GP correlation proposals
    gp_scale_smooth_U <- lapply(1:nclust, function(x) diag(2))

    # Get initial values--------------------------------------------
    alpha_cur <- init$alpha
    clust_probs_cur <- calc_class_probs(X_class, alpha_cur)
    clust_labels_cur <- init$clust_labels
    sigma_sq_cur <- init$sigma_sq
    B_cur <- init$B
    gp_scale_tilde_cur <- init$gp_scale_tilde
    gp_smooth_tilde_cur <- init$gp_smooth_tilde
    gp_scale_cur <- gp_scale_tilde_cur[clust_labels_cur]
    gp_smooth_cur <- gp_smooth_tilde_cur[clust_labels_cur]
    a_tilde_cur <- init$a_tilde
    b_tilde_cur <- init$b_tilde
    a_cur <- a_tilde_cur[clust_labels_cur]
    b_cur <- b_tilde_cur[clust_labels_cur]
    # If an initial value is not passed for clust_loc_tilde, don't include in model (set to 0 and don't sample)
    if(!is.null(init$clust_loc_tilde)){
      clust_loc_tilde_cur <- init$clust_loc_tilde
    }
    else{
      clust_loc_tilde_cur <- rep(0, nclust)
    }
    clust_loc_cur <- clust_loc_tilde_cur[clust_labels_cur]
    # If an initial value is passed for the beta_prior_var, then sample, otherwise fix
    beta_prior_var_cur <-ifelse(is.null(init$beta_prior_var),
                                     const$beta_prior_var, init$beta_prior_var)


    # Residual GP Inverse covariance--------------------------------------------
    gp_inv_cor_list_cur <- list()
    gp_cor_tilde_list_cur <- list()
    for(j in 1:nclust){
      cor_model_cur <-
        RandomFields::RMmatern(
          var = 1,
          scale = gp_scale_tilde_cur[j],
          nu = gp_smooth_tilde_cur[j],
          notinvnu = TRUE
        )
      gp_cor_tilde_list_cur[[j]] <-
        RandomFields::RFcovmatrix(cor_model_cur, x = coord)
    }
    for(j in 1:nrep){
      gp_inv_cor_list_cur[[j]] <-
        solveC(gp_cor_tilde_list_cur[[clust_labels_cur[j]]][subset_id_list[[j]],
                                                            subset_id_list[[j]]])
    }

    ##################################
    ######## Begin MCMC loop #########
    ##################################
    smp_loc <- 1
    if (progress_bar) {
      pb <- txtProgressBar(min = 0,
                           max = nsamples,
                           style = 3)
    }
    for (i in 1:niter) {
      if (progress_bar) {
        setTxtProgressBar(pb, smp_loc)
      }
      gamma1 <- c0 / (i + tune_k) ^ (c1)            # Automatic tuning constant
      pos <- (i - 1) %% win_len + 1
      nxt_pos <- i %% win_len + 1

      # Update sigma_sq --------------------------------------------------------
        sigma_sq_cur <-
          r_fc_sigma_sq_intensity(Y_list,
                                  clust_loc_cur,
                                  X_list,
                                  B_cur,
                                  gp_inv_cor_list_cur,
                                  a = a_cur,
                                  b = b_cur)

      # a, b -------------------------------------------------------------------
      a_tilde_cur <- update_a(
        a_tilde_cur,
        b_tilde_cur,
        sigma_sq_cur,
        unq_clust = 1:nclust,
        clust_labels = clust_labels_cur,
        lower_bound = const$a_lower_bound,
        upper_bound = const$a_upper_bound,
        b_prior_shape = const$b_prior_shape,
        b_prior_scale = const$b_prior_scale
      )
      a_cur <- a_tilde_cur[clust_labels_cur]

      b_tilde_cur <-
        update_b(
          b_tilde_cur,
          a_tilde_cur,
          sigma_sq_cur,
          unq_clust = 1:nclust,
          clust_labels = clust_labels_cur,
          prior_shape = const$b_prior_shape,
          prior_scale = const$b_prior_scale
        )
      b_cur <- b_tilde_cur[clust_labels_cur]

      if(!is.null(init$clust_loc_tilde)){
      # Update cluster location -------------------------------------------------
        clust_loc_tilde_cur <-
          r_fc_clust_loc_tilde(
            clust_loc_tilde_cur,
            Y_list,
            B_cur,
            X_list,
            gp_inv_cor_list_cur,
            sigma_sq_cur,
            clust_labels_cur,
            const$clust_loc_prior_var,
            const$prior_clust_loc_tilde_mean
          )
        clust_loc_cur <- clust_loc_tilde_cur[clust_labels_cur]
      }
      # Update beta ------------------------------------------------------------
      B_cur <-
        r_fc_beta_intensity(
          B_cur,
          D,
          Y_list,
          X_list,
          clust_loc_cur,
          gp_inv_cor_list_cur,
          sigma_sq_cur,
          prior_var = beta_prior_var_cur,
          use_recursive_prior = use_recursive_prior
        )

      # Update beta prior variance if init is passed ----------------------------
      if(!is.null(init$beta_prior_var)){
        beta_prior_var_cur <- r_fc_beta_prior_var(beta_prior_var = beta_prior_var_cur,
                                                  B_cur = B_cur,
                                                  D = D,
                                                  use_recursive_prior = use_recursive_prior)
      }

      # Update scale and smoothness parameters ---------------------------------
      if (i > thin_int * 5) {
        gp_scale_smooth_U <- list()
        for (k in 1:nclust) {
          gp_scale_smooth_cov <- cov(cbind(
            log(samples$gp_scale_tilde[,k]),
            log(samples$gp_smooth_tilde[,k])
          ),
          use = 'complete.obs')
          # Should try to come up with a better way to handle this
          if(all(gp_scale_smooth_cov ==0)){
            gp_scale_smooth_cov <- diag(2)*0.02
          }
          gp_scale_smooth_U[[k]] <- tryCatch(chol(gp_scale_smooth_cov),
                                             error = function(e) chol(diag(2)*0.02))
        }
      }

      gp_scale_smooth_tilde_updates <-
        update_scale_smooth_intensity(
          gp_scale_tilde_cur,
          prop_var$gp_scale_smooth_tilde[pos, ],
          Y_list,
          clust_loc_cur,
          X_list,
          B_cur,
          sigma_sq_cur,
          gp_smooth_tilde_cur,
          gp_inv_cor_list_cur,
          coord,
          subset_id_list,
          clust_labels_cur,
          1:nclust,
          const$gp_scale_range,
          const$gp_smooth_range,
          gp_scale_smooth_U
        )
      acpt$gp_scale_smooth_tilde[pos, ] <- gp_scale_smooth_tilde_updates$acpt
      for(k in 1:nclust){
        prop_var$gp_scale_smooth_tilde[nxt_pos, k] <-
          update_var(
            prop_var$gp_scale_smooth_tilde[pos, k],
            acpt_rt = mean(acpt$gp_scale_smooth_tilde[,k], na.rm = TRUE),
            opt_rt = opt_rt,
            gamma1 = gamma1
          )
      }
      gp_scale_tilde_cur <- gp_scale_smooth_tilde_updates$scale_tilde
      gp_smooth_tilde_cur <- gp_scale_smooth_tilde_updates$smooth_tilde
      gp_scale_cur <- gp_scale_tilde_cur[clust_labels_cur]
      gp_smooth_cur <- gp_smooth_tilde_cur[clust_labels_cur]
      gp_inv_cor_list_cur <- gp_scale_smooth_tilde_updates$gp_inv_cor_list
      for(j in 1:nclust){
        cor_model_cur <-
          RandomFields::RMmatern(
            var = 1,
            scale = gp_scale_tilde_cur[j],
            nu = gp_smooth_tilde_cur[j],
            notinvnu = TRUE
          )
        gp_cor_tilde_list_cur[[j]] <-
          RandomFields::RFcovmatrix(cor_model_cur, x = coord)
      }

      if(nclust > 1){   # When nclust is 1, the model is just a student T
        # Update alphas -----------------------------------------------------------
        alpha_updates <- update_alpha(alpha_cur, X_class, clust_labels_cur,
                                      prop_var$alpha[pos,,], const$alpha_prior_var)
        alpha_cur <- alpha_updates$A
        acpt$alpha[pos,,] <- alpha_updates$acpt
        clust_probs_cur <- calc_class_probs(X_class, alpha_cur)

        for(p in 1:p_class_covariates){
          for(k in 1:(nclust -1)){
            prop_var$alpha[nxt_pos, p, k] <-
              update_var(
                prop_var$alpha[pos, p, k],
                acpt_rt = mean(acpt$alpha[, p, k], na.rm = TRUE),
                opt_rt = opt_rt,
                gamma1 = gamma1
              )
          }
        }
        # Update cluster labels ---------------------------------------------------
        clust_labels_updates <- update_clust_labels(
          clust_labels_cur,
          unq_clust,
          clust_probs_cur,
          Y_list,
          clust_loc_tilde_cur,
          sigma_sq_cur,
          gp_inv_cor_list_cur,
          a_tilde_cur,
          b_tilde_cur,
          X_list,
          B_cur,
          gp_cor_tilde_list_cur,
          subset_id_list
        )
        clust_labels_cur <- clust_labels_updates$clust_labels_cur
        acpt$clust_labels[pos,] <- clust_labels_updates$acpt
        gp_inv_cor_list_cur <- clust_labels_updates$gp_inv_cor_list
        clust_loc_cur <- clust_loc_tilde_cur[clust_labels_cur]
        a_cur <- a_tilde_cur[clust_labels_cur]
        b_cur <- b_tilde_cur[clust_labels_cur]
        gp_scale_cur <- gp_scale_tilde_cur[clust_labels_cur]
        gp_smooth_cur <- gp_smooth_tilde_cur[clust_labels_cur]
      }

      # Save samples -----------------------------------------------------------
      if (i %% thin_int == 0) {
        samples$B[smp_loc, , ] <- B_cur
        samples$clust_loc_tilde[smp_loc, ] <- clust_loc_tilde_cur
        samples$sigma_sq[smp_loc, ] <- sigma_sq_cur
        samples$gp_scale_tilde[smp_loc,] <- gp_scale_tilde_cur
        samples$gp_smooth_tilde[smp_loc,] <- gp_smooth_tilde_cur
        samples$a_tilde[smp_loc,] <- a_tilde_cur
        samples$b_tilde[smp_loc,] <- b_tilde_cur
        samples$alpha[smp_loc,,] <- alpha_cur
        samples$clust_labels[smp_loc,] <- clust_labels_cur
        samples$beta_prior_var[smp_loc] <-  beta_prior_var_cur
        smp_loc <- smp_loc  + 1
      }
    }
    acpt_rts  <- list(
      gp_scale_tilde = colMeans(acpt$gp_scale_smooth_tilde, na.rm = TRUE),
      sigma_sq = colMeans(acpt$sigma_sq, na.rm = TRUE ),
      alpha = apply(acpt$alpha, c(2,3), mean, na.rm = TRUE),
      clust_labels = colMeans(acpt$clust_labels, na.rm = TRUE)
    )
    return(
      list(
        samples = samples,
        prop_var = prop_var,
        acpt_rts = acpt_rts,
        const = const
      )
    )
  }
