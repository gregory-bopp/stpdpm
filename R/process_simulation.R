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
rberngp <- function(nrep,
                    gvar,
                    gscl,
                    nu = 1 / 2,
                    mu = 0,
                    nloc = NULL,
                    locs = NULL) {
  if (is.null(locs)) {
    if (is.null(nloc)) {
      stop(
        "If an observation location matrix (locs) is not specified, the
        number of observation locations (nloc) must be passed."
      )
    }
    locs <- matrix(runif(2 * nloc), ncol = 2)
  }
  RandomFields::RFoptions(spConform = FALSE)
  model <-
    RandomFields::RMmatern(
      var = gvar,
      scale = gscl,
      nu = nu,
      notinvnu = TRUE
    )
  bin_gp <-
    RandomFields::RFsimulate(n = nrep,
                             model = model,
                             x = locs) + mu
  bin <- (bin_gp > 0)*1L
  # p <- 1/(1+ exp(-bin_gp))
  # bin <- p
  # bin[] <- rbinom(nloc*nrep, 1, p)
  return(list(
    bin = bin,
    bin_gp = bin_gp,
    locs = locs
  ))
}

#' Simulate T process
#' @description main function for simulating process.
#' @param nrep Number of replicate processes to simulate
#' @param gscl Matern covariance function scale parameter
#' @param nu Matern smoothness parameter
#' @param nloc (optional) number of locations to simulate. Must be passed if
#' \code{coord} matrix is not passed.
#' @param coord (optional) nloc x 2 matrix of spatial locations at which to
#' simulate process
#' @param mu mean of Gaussian process. Should be either a (1) scalar, (2) vector
#' of length nloc, or (3) matrix of dimension nloc x nrep
#' @param a T d.f.
#' @param b T scale
#' @return list
#' @export
#'
rstp <-   function(nrep,
                   a_tilde,
                   b_tilde,
                   gp_scale_tilde,
                   gp_smooth_tilde,
                   clust_loc_tilde = 0,
                   X_list = NULL,
                   B = NULL,
                   nloc = NULL,
                   coord = NULL,
                   clust_labels = NULL,
                   X_class = NULL,
                   alpha = NULL,
                   subset_id_list = NULL,
                   pred_coord = NULL,
                   X_pred = NULL) {
  if (is.null(coord)) {
    if (is.null(nloc)) {
      stop(
        "If an observation location matrix (coord) is not specified, the
        number of observation locations (nloc) must be passed."
      )
    }
    coord <- matrix(runif(2 * nloc), ncol = 2)
  }
  else{
    nloc <- nrow(coord)
  }
  if(is.null(subset_id_list)){
    subset_id_list <- lapply(1:nrep, function(x) 1:nloc)
  }
  else{
    if(!all(sapply(subset_id_list, function(x) all(sort(x) ==x)))){
      stop("Elements of subset_id_list should be in increasing order with
           no duplicates.")
    }
  }
  # Check for prediction grid
  if(!is.null(pred_coord)){
    all_coord <- rbind(coord, as.matrix(pred_coord))
    nloc_pred <- nrow(pred_coord)
    nloc_all <- nloc + nloc_pred
  }
  else{
    all_coord <- coord
    nloc_all <- nloc
  }


  if (is.null(clust_labels)) {
    if ((is.null(alpha)) || (is.null(X_class))) {
      stop("Must pass one of clust_labels, or X_class and alpha")
    }
    probs <- calc_class_probs(X_class, alpha)
    nclust <- ncol(probs)
    clust_labels <- apply(probs, 1, function(x) {
      sample(1:nclust,
             size = 1,
             replace = F,
             prob = x)
    })
  }
  else{
    nclust <- length(unique(clust_labels))
  }
  expand_scalars(
    c(
      "a_tilde",
      "b_tilde",
      "gp_scale_tilde",
      "gp_smooth_tilde",
      "clust_loc_tilde"
    ),
    nclust
  )
  a <- a_tilde[clust_labels]
  b <- b_tilde[clust_labels]
  gp_scale <- gp_scale_tilde[clust_labels]
  gp_smooth <- gp_smooth_tilde[clust_labels]
  clust_loc <- clust_loc_tilde[clust_labels]
  # Simulating an inverse gamma with shape, scale parameters a, b is equivalent
  # to taking the reciprocal of a gamma with shape, scale parameters a, 1/b
  sigma_sq <- 1 / rgamma(nrep, shape = a / 2, rate = a * b / 2)
  RandomFields::RFoptions(spConform = FALSE)
  gp <- matrix(NA_real_, nrow = nloc_all, ncol = nrep)
  for (i in 1:nrep) {
    model <-
      RandomFields::RMmatern(
        var = 1,
        scale = gp_scale[i],
        nu = gp_smooth[i],
        notinvnu = TRUE
      )
    gp[, i] <-  RandomFields::RFsimulate(n = 1,
                                   model = model,
                                   x = all_coord)
  }
  Y_all <- gp %*% diag(sqrt(sigma_sq))
  Y <- Y_all[1:nloc,]
  Y_list <- list()
  if(is.null(X_list)){
    for(i in 1:nrep){
      Y_list[[i]] <- Y[subset_id_list[[i]],i] + clust_loc[i]
    }
  }
  else{
    for(i in 1:nrep){
      Y_list[[i]] <- Y[subset_id_list[[i]],i] + clust_loc[i] + X_list[[i]] %*% B[,i]
    }
  }
  # Holdout grid
  if(!is.null(pred_coord)){
    if(is.null(X_pred)){
      Y_hold <- Y_all[(nloc+1):(nloc_all),] + +
        matrix(clust_loc, nrow = nloc_pred, ncol = nrep, byrow = T)
    }
    else{
      Y_hold <- Y_all[(nloc+1):(nloc_all),] +
        matrix(clust_loc, nrow = nloc_pred, ncol = nrep, byrow = T) +
        X_pred %*% B
    }
  }
  else{
    Y_hold <- NULL
  }

  return(
    list(
      Y_list = Y_list,
      Y_hold = Y_hold,
      B = B,
      a = a,
      b = b,
      gp_scale = gp_scale,
      gp_smooth = gp_smooth,
      sigma_sq = sigma_sq,
      gp = gp,
      coord = coord,
      clust_labels = clust_labels,
      a_tilde = a_tilde,
      b_tilde = b_tilde,
      gp_scale_tilde = gp_scale_tilde,
      gp_smooth_tilde = gp_smooth_tilde,
      clust_loc_tilde = clust_loc_tilde,
      alpha = alpha,
      subset_id_list = subset_id_list
    )
  )
}


