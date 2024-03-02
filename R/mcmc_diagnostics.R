#' @export
intensity_traceplots <-
  function(samples,
           which_reps = c(1, dim(samples$B)[3]),
           fig_filepath = "./intensity_trace.pdf") {
    pdf(fig_filepath)
    nrep <- dim(samples$B)[3]
    nclust <- ncol(samples$gp_scale_tilde)
    nmcmc <- dim(samples$B)[1]

    #  Counts of cluster labels
    props <- matrix(NA_real_, nmcmc, nclust)
    for(k in 1:nclust){
      props[,k] <- apply(out$samples$clust_labels, 1, function(x) sum(x == k)/length(x))
    }
    par(mfrow =c(1,1))
    matplot(props, type = "l", lty = 1,
            xlab = "iteration", ylab = "number of replicates in cluster",
            ylim = c(0,max(props)))
    abline(h = 0, lty = 2)
    legend("topright", legend = 1:nclust, col = 1:nclust, lty = 1)

    props <- matrix(NA_real_, nrep, nclust)
    for(k in 1:nclust){
      props[,k] <- apply(out$samples$clust_labels, 2, function(x) sum(x == k)/length(x))
    }
    matplot(props, lty = 1, pch = 21,
            xlab = "time", ylab = "prop samples of cluster k")
    legend("topright", legend = 1:nclust, col = 1:nclust, lty = 1)

    # Beta prior variance
    plot(samples$beta_prior_var,
          type = "l",
           main = "Beta prior variance",
           ylab = "sample")

    # Betas
    par(mfrow = c(2, 2))
    for (p in 1:dim(samples$B)[2]) {
      for (i in which_reps) {
        plot(
          samples$B[, p, i],
          type = "l",
          main = paste0("beta_", p, " time = ", i),
          ylab = "sample"
        )
      }
    }

    # cluster location parameters
    for (k in 1:nclust) {
      plot(
        samples$clust_loc_tilde[, k],
        type = "l",
        main = paste0("clust_loc_tilde_", k),
        ylab = "sample"
      )
    }

    # gp scale tilde
    for (k in 1:nclust) {
      plot(
        samples$gp_scale_tilde[, k],
        type = "l",
        main = paste0("gp_scale_tilde_", k),
        ylab = "sample"
      )
    }

    # gp smooth tilde
    for (k in 1:nclust) {
      plot(
        samples$gp_smooth_tilde[, k],
        type = "l",
        main = paste0("gp_smooth_tilde_", k),
        ylab = "sample"
      )
    }

    # a tilde
    for (k in 1:nclust) {
      plot(
        samples$a_tilde[, k],
        type = "l",
        main = paste0("a_tilde_", k),
        ylab = "sample"
      )
    }

    # b tilde
    for (k in 1:nclust) {
      plot(
        samples$b_tilde[, k],
        type = "l",
        main = paste0("b_tilde_", k),
        ylab = "sample"
      )
    }

    # alpha
    if(nclust!=1){
      for (p in 1:dim(samples$alpha)[2]) {
        for (k in 1:(nclust - 1)) {
          plot(
            samples$alpha[, p, k],
            type = "l",
            main = paste0("alpha_", p, " cluster: ", k),
            ylab = "sample"
          )
        }
      }
    }

    # sigma_sq
    for (i in which_reps) {
      plot(
        samples$sigma_sq[, i],
        type = "l",
        main = paste0("sigma_sq_", i),
        ylab = "sample"
      )
    }

    # cluster labels
    if(nclust!=1){
      for (i in which_reps) {
        plot(
          samples$clust_labels[, i],
          type = "l",
          main = paste0("cluster ", i),
          ylab = "sample"
        )
      }
    }
    dev.off()
  }

#' @export
intensity_truth_corplot <-
  function(samples, truth, fig_filepath = "./intensity_truth_cor.pdf") {
    pdf(fig_filepath)
    par(mfrow = c(2, 2))
    nclust <- length(unique(truth$clust_labels))
    if(nclust > 1){
      # Cluster labels
      clmn <- apply(samples$clust_labels, 2, mean)
      plot(truth$clust_labels,
           clmn , main = "")
      title(main = paste("cor: ", round(cor(
        truth$clust_labels, clmn
      ), 2)))
      abline(a = 0, b = 1)

      # alpha
      alphamn <- apply(samples$alpha, c(2,3), mean)
      plot(c(truth$alpha), c(alphamn))
      title(main = paste("cor: ", round(cor(
        c(truth$alpha), c(alphamn)
      ), 2)))
      abline(a = 0, b = 1)
    }
    # gp_scale_tilde and gp_smooth_tilde
    gp_scale_tilde_pm <- apply(samples$gp_scale_tilde, 2, mean)
    gp_smooth_tilde_pm <- apply(samples$gp_smooth_tilde, 2, mean)
    plot(truth$gp_scale_tilde, gp_scale_tilde_pm )
    abline(a = 0, b = 1)
    title(main = paste("cor: ", round(cor(
      gp_scale_tilde_pm, truth$gp_scale_tilde
    ), 2)))
    plot(truth$gp_smooth_tilde, gp_smooth_tilde_pm )
    abline(a = 0, b = 1)
    title(main = paste("cor: ", round(cor(
      gp_smooth_tilde_pm, truth$gp_smooth_tilde
    ), 2)))
    # a and b
    apm <- apply(samples$a_tilde, 2, mean)
    bpm <- apply(samples$b_tilde, 2, mean)
    plot(truth$a_tilde, apm)
    abline(a = 0, b = 1)
    title(main = paste("cor: ", round(cor(
      apm, truth$a_tilde
    ), 2)))
    plot(truth$b_tilde, bpm)
    abline(a = 0, b = 1)
    title(main = paste("cor: ", cor(bpm, truth$b_tilde)))

    # cluster location
    clocmn <- apply(samples$clust_loc_tilde, 2, mean)
    plot(truth$clust_loc_tilde,
         clocmn)
    abline(a = 0, b = 1)
    title(main = paste("cor: ", round(cor(
      clocmn, truth$clust_loc_tilde
    ), 2)))

    # beta
    for (p in 1:nrow(truth$B)) {
      Bmn <- apply(samples$B, c(2, 3), mean)
      plot(
        truth$B[p, ],
        Bmn[p, ],
        xlab = paste0("true beta_", p),
        ylab = paste0("Post. Mean beta_", p)
      )
      abline(a = 0, b = 1)
      title(main = paste("cor: ", round(cor(Bmn[p, ], B[p, ]), 2)))
    }

    # sigma_sq
    ssmn <- apply(samples$sigma_sq, 2, mean)
    plot(truth$sigma_sq, ssmn)
    abline(a = 0, b = 1)
    title(main = paste("cor: ",  round(cor(
      ssmn, truth$sigma_sq
    ), 2)))

    dev.off()
  }


#' Traceplots for occurrence model
#'
#' @param samples   mcmc list of samples
#' @param which_locs which locations to plot for latent Z samples
#' @param which_reps which time replicates to plot
#' @param fig_filepath filepath to store traceplots
#' @export
occur_traceplots <- function(samples,
                             which_locs = 1,
                             which_reps = c(1, dim(samples$beta)[3]),
                             fig_filepath = "./occur_trace.pdf") {
  pdf(fig_filepath)
  nrep <- dim(samples$beta)[3]
  par(mfrow = c(2, 2))
  # Correlation parameters
  plot(samples$gp_scale,
       type = "l",
       main = "gp_scale",
       ylab = "sample")
  plot(samples$gp_smooth,
       type = "l",
       main = "gp_smooth",
       ylab = "sample")

  # Beta prior variance
  plot(samples$beta_prior_var,
       type = "l",
       main = "Beta prior variance",
       ylab = "sample")

  # Betas
  for (p in 1:dim(samples$beta)[2]) {
    for (i in which_reps) {
      plot(
        samples$beta[, p, i],
        type = "l",
        main = paste0("beta_", p, " time = ", i),
        ylab = "sample"
      )
    }
  }
  # Z
  for (i in which_reps) {
    for (j in which_locs) {
      plot(
        samples$Z[, j, i],
        type = "l",
        main = paste0("Z time:", i, " loc:", j),
        ylab = "sample"
      )
    }
  }
  dev.off()
}


#' @export
occur_truth_corplot <-
  function(samples, truth, fig_filepath = "./occur_truth_cor.pdf") {
    pdf(fig_filepath)
    par(mfrow = c(2, 2))
    # Correlation parameters
    plot(samples$gp_scale,
         type = "l",
         main = "gp_scale",
         ylab = "sample")
    abline(h = truth$gp_scale)
    plot(samples$gp_smooth,
         type = "l",
         main = "gp_smooth",
         ylab = "sample")
    abline(h = truth$gp_smooth)

    # Betas
    bmn <- apply(samples$beta, c(2, 3), mean)
    plot(truth$beta, bmn)
    abline(a = 0, b = 1)
    title(main = paste("cor: ", round(cor(c(bmn), c(truth$beta)), 2)))

    # Latent Z
    Zmn <- apply(samples$Z, c(2, 3), mean)
    plot(truth$Z, Zmn)
    abline(a = 0, b = 1)
    title(main = paste("cor: ", round(cor(c(Zmn), c(truth$Z)), 2)))

    dev.off()
  }



#' @export
ppred_intensity_quantiles_one <- function(quantiles,
                                          which_reps,
                                          nclust,
                                          subset_id_list,
                                          obs_coords,
                                          pred_coords,
                                          Y_list,
                                          X_list_obs,
                                          X_pred,
                                          B,
                                          clust_loc,
                                          sigma_sq,
                                          gp_scale,
                                          gp_smooth,
                                          clust_labels,
                                          flatten = TRUE) {
  Y_quantiles <-
    array(NA_real_, dim = c(nrow(pred_coords), length(which_reps), length(quantiles)))
  models <- list()
  for (k in 1:nclust) {
    models[[k]] <- RandomFields::RMmatern(
      var = 1,
      scale = gp_scale[k],
      nu = gp_smooth[k],
      notinvnu = TRUE
    )
  }
  npred_coord <- nrow(pred_coords)
  Yhat_mean <-
    matrix(NA_real_, nrow = npred_coord, ncol = length(which_reps))
  l <- 0
  for (i in which_reps) {
    l <- l + 1
    eps <-
      (Y_list[[i]] - clust_loc[clust_labels[i]] - X_list_obs[[i]] %*% B[, i]) /
      sqrt(sigma_sq[i])
    eps_mean <-
      RFinterpolate(
        models[[clust_labels[i]]],
        x = pred_coords[, 1],
        y = pred_coords[, 2],
        data = data.frame(x = obs_coords[subset_id_list[[i]], 1],
                          y = obs_coords[subset_id_list[[i]], 2], G = eps)
      )
    Yhat_mean[, l] <- eps_mean + clust_loc[clust_labels[i]] +
                        X_pred %*% B[, i]
    for (j in 1:length(quantiles)) {
      Y_quantiles[, l, j] <-
        qnorm(quantiles[j], Yhat_mean[, l], sd = sqrt(sigma_sq[i]))
    }
  }
  if (flatten) {
    return(c(Y_quantiles))
  }
  return(Y_quantiles)
}

#' @export
ppred_intensity_quantiles_for_trace <- function(samples,
                                                thin_int,
                                                quantiles,
                                                which_reps,
                                                nclust,
                                                subset_id_list,
                                                obs_coords,
                                                pred_coords,
                                                Y_list,
                                                X_list_obs,
                                                X_pred){
  nmcmc <- nrow(samples$clust_labels)
  mcmc_its <- seq(1, nmcmc, by = thin_int)
  i <- 0
  qhat <- matrix(NA_real_, nrow = length(mcmc_its),
                 ncol = length(quantiles)*length(which_reps)*nrow(pred_coords))
  for(mcmc_it in mcmc_its){
    i <- i + 1
    qhat[i,] <- ppred_intensity_quantiles_one(
      quantiles = quantiles,
      which_reps = which_reps,
      nclust = nclust,
      subset_id_list = subset_id_list,
      obs_coords = obs_coords,
      pred_coords = pred_coords,
      Y_list = Y_list,
      X_list_obs = X_list_obs,
      X_pred = X_pred,
      B = samples$B[mcmc_it, ,],
      clust_loc = samples$clust_loc_tilde[mcmc_it,],
      sigma_sq = samples$sigma_sq[mcmc_it, ],
      clust_labels = samples$clust_labels[mcmc_it,],
      gp_scale = samples$gp_scale_tilde[mcmc_it,],
      gp_smooth = samples$gp_smooth_tilde[mcmc_it,],
      flatten = TRUE
    )
  }
  return(qhat)
}
