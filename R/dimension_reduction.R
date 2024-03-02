
# n_near_neighbors <- 20              # Number of near-neighbors to use
# p_reduced_dim <- 10                 # Size of reduced dimension to use

#' Laplace Eigenmap Dimension Reduction (see Belkin, 2001)
#'
#' @param dists distance matrix that is used to create adjacency matrix. e.g.
#' the i,jth entry gives the distance between the covariates for the ith time
#' period and jth time period.
#' @param n_near_neighbors number of near neighbors to use when creating
#'                         adjacency matrix
#' @param p_reduced_dim reduced dimension of vectors. If NULL, will return full
#'            laplace eigenmap feature set. The first dimension of the Laplace
#'            eigenmap is a constant (eigenvalue = 0), and so is ommitted in the
#'            returned matrix unless p_reduced_dim = NULL
#' @return Laplace eigenmap embedding of observations. e.g. the ith time point
#'         in the embedding space will be represented as the ith row of the
#'         output.
#' @export
lapeigenmap <- function(dists, n_near_neighbors, p_reduced_dim = NULL) {
  G <- t(apply(dists, 1,
               function(x) {
                 1L * (x %in%
                         sort(x,
                              partial = n_near_neighbors)[1:n_near_neighbors])
               }))                 # Adjacency matrix of nearest neighbors
  W <- pmax(G, t(G))               # Weight matrix needs to be symmetric
  D <- diag(rowSums(W))            # Diagonal weight matrix(Belkin & Niyogi, 2001)
  L <- D - W                       # Laplacian matrix
  lamap <- geigen::geigen(L, D,    # Generalized eigenvalue problem
                          symmetric = TRUE,
                          only.values = FALSE)
  if(is.null(p_reduced_dim)){
    return(lamap$vectors)
  }
  return(lamap$vectors[, 2:(1 + p_reduced_dim)])
}
