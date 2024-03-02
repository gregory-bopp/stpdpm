#' @export
gauss_kernel_gram <- function(X, theta = 1.0){
  return(gauss_kernel_gramC(X, theta))
}

#' @export
calc_dist <- function(X, Y = NULL){
  if(is.null(Y)){
    return(calc_dist_selfC(X))
  }
  else{
    return(calc_distC(X, Y))
  }
}
