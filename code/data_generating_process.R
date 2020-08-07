# Data-generating Process

library(mvtnorm)



generate_data <- function(n=200, n_F_attr, n_G_attr, n_H_attr, treatment_effect,
                          unconfoundedness_rate) {


  # Test function arguments
  fraction <- n_H_attr / 3
  if (fraction%%1 != 0) {
    print("Warning: Number of house attributes must be divisible by 3.")
  }



  # Draw features explaining treatment which are independent of outcome
  F <- replicate(n=n_F_attr, rnorm(n=n, mean=2, sd=1))

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





  # Data-generating process

  # Degree of sparseness in attributes explaining treatment only
  sparsity_rate_F <- .75
  beta_F <- rep(0, dim(F)[2])
  zero_coefs_F <- sample(1:dim(F)[2], size=sparsity_rate_F*dim(F)[2])
  beta_F[zero_coefs_F] <- 0


  # "Unconfoundedness rate" :
  # Sparsity in both models explaining treatment and outcome,
  # i.e. the same geographical features are important to predict D and Y.
  # However, different effect sizes are allowed.
  sparsity_rate_GDY <- unconfoundedness_rate
  zero_coefs_GDY <- sample(1:dim(G)[2], size=sparsity_rate_GDY*dim(G)[2])

  beta_GD <- rep(.25, dim(G)[2])
  beta_GD[zero_coefs_GDY] <- 0

  beta_GY <- rep(.25, dim(G)[2])
  beta_GY[zero_coefs_GDY] <- 0



  # Degree of sparseness in housing attributes explaining house prices
  sparsity_rate_H <- .75
  beta_H <- rep(1, dim(H)[2])
  zero_coefs_H <- sample(1:dim(H)[2], size=sparsity_rate_H*dim(H)[2])
  beta_H[zero_coefs_H] <- 0


  # Collect information which features' coefficients are set to zero
  sparsity_identifier_F <- rep(0, dim(F)[2])
  sparsity_identifier_G <- rep(0, dim(G)[2])
  sparsity_identifier_H <- rep(0, dim(H)[2])
  sparsity_identifier_F[zero_coefs_F] <- 1
  sparsity_identifier_G[zero_coefs_GDY] <- 1
  sparsity_identifier_H[zero_coefs_H] <- 1
  sparsity_identifier <- append(sparsity_identifier_F, sparsity_identifier_G)
  sparsity_identifier <- append(sparsity_identifier, sparsity_identifier_H)
  # Collect identifier for variables which are explaining either D or Y
  true_covariate_identifier <- seq(from=1, to=(n_F_attr+n_G_attr+n_H_attr), by=1)
  true_covariate_identifier[sparsity_identifier==0] <- NaN
  true_covariate_identifier <- seq(from=1, to=(n_G_attr+n_H_attr), by=1)
  true_covariate_identifier[sparsity_identifier==0] <- NaN


  # Treatment effect
  beta_DY <- treatment_effect

  # Add white noise to both models for treatment and outcome
  eps_D <- rnorm(n=n, mean=0, sd=1)
  eps_Y <- rnorm(n=n, mean=0, sd=1)

  # Treatment
  # --- linear model
  D <- 0 + G %*% beta_GD + eps_D
  D <- 0 + G %*% beta_GD + F %*% beta_F + eps_D
  # --- logistic model
  #proba_binom_model <- exp(0 + G %*% beta_GD + eps_D) / (1 + exp(0 + G %*% beta_GD + eps_D))
  #D <- cbind(rbinom(n=n, size=1, prob=proba_binom_model))

  # Outcome
  y <- 0 + D %*% beta_DY + G %*% beta_GY + H %*% beta_H + eps_Y



  # Create vector with covariate names.
  colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
  colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
  colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))

  # Required for first selection step
  colnames_dataset <- c("y", "D", colname_F, colname_G, colname_H)

  data <- data.frame(cbind(y, D, F, G, H))
  names(data) <- colnames_dataset
  #return(list(y=y, D=D, F=F, G=G, H=H, true_covariate_identifier=true_covariate_identifier))
  return(list(data=data, true_covariate_identifier=true_covariate_identifier))
}


















#
