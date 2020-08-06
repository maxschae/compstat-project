### Least Squares function

least_squares_estimator <- function(y, X) {
  # Computes the least squares estimator including the intercept.
  
  # Args:
  #     y (vector): dependent variable
  #     X (matrix): independent variables
  #
  # Return (list)
  #     beta: OLS estimator

  X <- cbind(rep(1, dim(X)[1]), X)
  k <- dim(X)[2]

  if (qr(X)$rank == k) {
    # Estimate coefficients
    beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
  }
  else if (qr(X)$rank != dim(X)[2]) {
    print("Warning: covariate matrix has not full column rank")
  }

  return(beta_hat)
}
