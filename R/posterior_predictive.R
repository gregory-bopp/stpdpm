#' Make posterior predictive draws from occurrence model at new locations
#'
#' @param obs_coords number of obs loc x 2 matrix of observation coordinates
#' @param pred_coords number of prediction loc x 2 matrix of prediction
#'   coordinates
#' @param scale Matern scale parameter
#' @param smooth Matern smoothness parameter
#' @param Z Conditional GP at observation coordinates (from mcmc)
#' @param ndraws NOT IMPLEMENTED YET. number of draws to make (default = 1). If ndraws > 1,
#' independent draws will be stored columnwise
#'
#' @return vector or matrix of posterior predictive draws of GP
# output dimension: [loc, time, replicate draw]
#' @export
ppred_occur_one_iter <-
  function(obs_coords,
           pred_coords,
           scale,
           smooth,
           Z,
           Xarr_obs,
           Xarr_pred,
           B,
           which_reps = NULL,
           ndraws = 1,
           return_binary = FALSE) {
    RandomFields::RFoptions(pch = "")
    nrep <- ncol(Z)
    nloc <- nrow(pred_coords)
    if(is.null(which_reps)){
      which_reps <- 1:nrep
    }
    model <- RandomFields::RMmatern(
      var = 1,
      scale = scale,
      nu = smooth,
      notinvnu = TRUE
    )
    RandomFields::RFoptions(spConform = FALSE)
    Zhat <- array(NA_real_, dim = c(nloc, length(which_reps), ndraws))
    l <- 0
    for(i in which_reps){
      l <- l + 1
      ep <- Z[,i] - Xarr_obs[,,i] %*%B[,i]
      Zhat[,l,] <- RandomFields::RFsimulate(
        n = ndraws,
        model = model,
        x = pred_coords[, 1],
        y = pred_coords[, 2],
        data = data.frame(x = obs_coords[, 1], y = obs_coords[, 2], G = ep)
      ) + c(Xarr_pred[,,i] %*% B[,i])
    }
    if(return_binary){
      return(1L*(Zhat > 0))
    }
    return(Zhat)
  }


#' Make posterior predictive draws from occurrence model
#' @param samples samples object from mcmc_occur
#' @param ndraws number of replicate draws to make per time and mcmc iteration
#' @param which_reps which time replicates to samples
#' @param obs_coords n1x2 matrix of observation coordinates (should match scale model was fit on)
#' @param pred_coords n2x2 matrix of prediction grid coordinates (should match scale model was fit on)
#' @param Xarr_obs [nloc, nvar, nrep] array of covariates for observation locations
#' @param Xarr_pred [nloc, nvar, nrep] array of covariates for prediction locations
#' @param thin_int thinning interval (if not all mcmc iterations should be used. Don't also pass mcmc_its)
#' @param mcmc_its vector of mcmc_iterations if not all should be used
#' @param return_binary Should just the binary process be returned, or should the full latent GP be returned?
#' @return array of posterior predictive draws of binary process
# output dimension: [mcmc, loc, time, replicate draw]
#' @export
ppred_occur <- function(samples,
                        ndraws,
                        which_reps,
                        obs_coords,
                        pred_coords,
                        Xarr_obs,
                        Xarr_pred,
                        thin_int = 1,
                        mcmc_its = NULL,
                        return_binary = FALSE){
  if(is.null(mcmc_its)){
    nmcmc <- nrow(samples$beta)
    mcmc_its <- seq(1, nmcmc, by = thin_int)
    nmcmc <- length(mcmc_its)
  }
  else{
    nmcmc <- length(mcmc_its)
  }
  i <- 0
  yhat <- array(NA_real_, dim = c(nmcmc, nrow(pred_coords), length(which_reps), ndraws))
  for(mcmc_it in mcmc_its){
    i <- i + 1
    yhat[i, , , ] <-
      ppred_occur_one_iter(
        obs_coords = obs_coords,
        pred_coords = pred_coords,
        scale = samples$gp_scale[i],
        smooth = samples$gp_smooth[i],
        Z = samples$Z[i,,],
        Xarr_obs = Xarr_obs,
        Xarr_pred = Xarr_pred,
        B = samples$beta[i,,],
        which_reps = which_reps,
        ndraws = ndraws,
        return_binary = return_binary
      )
  }
  return(yhat)
}


# Output format: [loc, time, replicate draw]
#' @export
ppred_intensity_one_iter <-function(ndraws,
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
                          clust_labels) {
  RandomFields::RFoptions(pch = "")
  Yhat <-
    array(NA_real_, dim = c(nrow(pred_coords), length(which_reps), ndraws))
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
  l <- 0
  for (i in which_reps) {
    l <- l + 1
    eps <-
      (Y_list[[i]] - clust_loc[clust_labels[i]] - X_list_obs[[i]] %*% B[, i]) /
        sqrt(sigma_sq[i])
    epshat <-
      sqrt(sigma_sq[i])*RandomFields::RFsimulate(
        n = ndraws,
        model = models[[clust_labels[i]]],
        x = pred_coords[, 1],
        y = pred_coords[, 2],
        data = data.frame(x = obs_coords[subset_id_list[[i]], 1],
                          y = obs_coords[subset_id_list[[i]], 2], G = eps)
      )
    Yhat[, l, ] <- epshat + matrix(clust_loc[clust_labels[i]] +
                                X_pred %*% B[, i],
                                nrow = npred_coord, ncol = ndraws)
  }
  return(Yhat)
}


# Output format: [mcmc, loc, time, replicate draw]
#' @export
ppred_intensity <- function(samples,
                            ndraws,
                            which_reps,
                            nclust,
                            subset_id_list,
                            obs_coords,
                            pred_coords,
                            Y_list,
                            X_list_obs,
                            X_pred,
                            thin_int = 1,
                            mcmc_its = NULL){
  if(is.null(mcmc_its)){
    nmcmc <- nrow(samples$clust_labels)
    mcmc_its <- seq(1, nmcmc, by = thin_int)
    nmcmc <- length(mcmc_its)
  }
  else{
    nmcmc <- length(mcmc_its)
  }
  i <- 0
  yhat <- array(NA_real_, dim = c(nmcmc, nrow(pred_coords), length(which_reps), ndraws))
  for(mcmc_it in mcmc_its){
    i <- i + 1
    yhat[i,,,] <- ppred_intensity_one_iter(
      ndraws,
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
      gp_smooth = samples$gp_smooth_tilde[mcmc_it,]
    )
  }
  return(yhat)
}
