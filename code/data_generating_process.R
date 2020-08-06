# Data-generating Process

library(mvtnorm)



generate_data <- function(n=200, n_F_attr, n_G_attr, n_H_attr) {


  # Test function arguments
  fraction <- n_H_attr / 3
  if (fraction%%1 != 0) {
    print("Warning: Number of house attributes must be divisible by 3.")
  }

  # Draw house attributes, continuous and integers
  # --- correlated
  n_H_attr_cont <- n_H_attr / 3
  n_H_attr_int <- n_H_attr / 3
  n_H_attr_cat <- n_H_attr / 3

  H_mean <- rep(0, n_H_attr_cont)

  A <- matrix(rnorm(n=n_H_attr_cont**2, mean=0, sd=1), ncol=n_H_attr_cont)  # draw random matrix
  H_sigma <- t(A) %*% A  # ensure matrix is positive-semidefinite
  #TODO sigma matrix must have positive covariances, see code snippets
  # also: ensure more correlation

  #TODO, for now, Hs are uncorrelated --- dgp_1
  H_variances <- runif(n=n_H_attr_cont, min=1, max=4)
  H_sigma <- diag(H_variances)

  H_continuous <- rmvnorm(n=n, mean=H_mean, sigma=H_sigma)
  H_integers <- floor(rmvnorm(n=n, mean=H_mean, sigma=H_sigma))

  # or take subset of H_integers and one-hot encode them
  H_categorical <- replicate(n=n_H_attr_cat, rbinom(n=n, size=3, prob=.5))
  #TODO one-hot encoding

  H <- cbind(H_continuous, H_integers, H_categorical)



  # Draw geographical features
  # --- uncorrelated
  G_mean <- rep(0, n_G_attr/2)

  G_variances <- runif(n=n_G_attr/2, min=1, max=4)
  G_sigma <- diag(G_variances)

  G_continuous <- rmvnorm(n=n, mean=G_mean, sigma=G_sigma)
  G_categorical <- floor(rmvnorm(n=n, mean=G_mean, sigma=G_sigma))
  #TODO one-hot encode subset (or all) of G_categorical vars

  # Using a binomial distribution.
  #G_categorical <- replicate(n=n_G_attr/2, rbinom(n=n, size=10, prob=.5))

  G <- cbind(G_continuous, G_categorical)


  # Draw features explaining treatment which are independent of outcome
  F <- replicate(n=n_F_attr, rnorm(n=n, mean=2, sd=1))


  # Data-generating process
  eps_D <- rnorm(n=n, mean=0, sd=1)
  eps_Y <- rnorm(n=n, mean=0, sd=1)

  #TODO
  # Degree of sparseness in geographical features explaining treatment
  #sparsity_rate_GD <- .5
  #beta_GD <- rep(1.5, dim(G)[2])
  #zero_coefs <- sample(1:dim(G)[2], size=sparsity_rate_GD*dim(G)[2])
  #beta_GD[zero_coefs] <- 0

  # Degree of sparseness in geographical features explaining house prices
  #TODO Note that beta_GY << beta_GD
  #sparsity_rate_GY <- .5
  #beta_GY <- rep(.15, dim(G)[2])
  #zero_coefs <- sample(1:dim(G)[2], size=sparsity_rate_GY*dim(G)[2])
  #beta_GY[zero_coefs] <- 0


  # "Unconfoundedness rate" :
  # Sparsity in both models explaining treatment and outcome,
  # i.e. the same geographical features are important to predict D and Y.
  # However, different effect sizes are allowed.
  sparsity_rate_GDY <- .5
  zero_coefs_GDY <- sample(1:dim(G)[2], size=sparsity_rate_GDY*dim(G)[2])

  beta_GD <- rep(1, dim(G)[2])
  beta_GD[zero_coefs_GDY] <- 0

  beta_GY <- rep(1, dim(G)[2])
  beta_GY[zero_coefs_GDY] <- 0



  # Degree of sparseness in housing attributes explaining house prices
  sparsity_rate_H <- .5
  beta_H <- rep(1, dim(H)[2])
  zero_coefs <- sample(1:dim(H)[2], size=sparsity_rate_H*dim(H)[2])
  beta_H[zero_coefs] <- 0


  # Degree of sparseness in further (geographical) attributes explaining treatment
  sparsity_rate_F <- .5
  beta_F <- rep(.5, dim(F)[2])
  zero_coefs <- sample(1:dim(F)[2], size=sparsity_rate_H*dim(F)[2])
  beta_F[zero_coefs] <- 0


  # Treatment effect
  beta_DY <- -1


  # Treatment
  # --- linear model
  D <- 0 + G %*% beta_GD + eps_D
  D <- 0 + G %*% beta_GD + F %*% beta_F + eps_D
  # --- logistic model
  #proba_binom_model <- exp(0 + G %*% beta_GD + eps_D) / (1 + exp(0 + G %*% beta_GD + eps_D))
  #D <- cbind(rbinom(n=n, size=1, prob=proba_binom_model))

  # Outcome
  y <- 0 + D %*% beta_DY + G %*% beta_GY + H %*% beta_H + eps_Y


  return(list(y=y, D=D, F=F, G=G, H=H))
}


















#
