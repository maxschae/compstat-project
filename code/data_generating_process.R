### CompStat project --- Max Schäfer
# Define the data-generating processes


generate_data <- function(n=200, n_F_attr=30, n_G_attr=30, n_H_attr=30, treatment_effect=.5,
                          beta_GD_size=.5, beta_GY_size=.25, beta_H_size=1, beta_F_size=1,
                          unconfoundedness_rate=.9) {


  # Test function arguments
  fraction <- n_H_attr / 3
  if (fraction%%1 != 0) {
    print("Warning: Number of house attributes must be divisible by 3.")
  }




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

  #TODO
  G <- replicate(n=n_G_attr, rnorm(n=n, mean=0, sd=1))


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



  #TODO Simplest case with i.i.d N(0, 1) variables and beta=1
  #TODO
  #TODO HEADS UP WITH F --- Collider bias incoming?
  F <- replicate(n=n_F_attr, rnorm(n=n, mean=0, sd=1))
  G <- replicate(n=n_G_attr, rnorm(n=n, mean=0, sd=1))
  H <- replicate(n=n_H_attr, rnorm(n=n, mean=0, sd=1))


  # Draw positively correlated regressors
  sigma_G <- matrix(.01, nrow=n_G_attr, ncol=n_G_attr)

  for (j in 1:(n_G_attr)) {
    power <- 1 - j
    a <- seq(from=power, to=(n_G_attr-j), by=1)
    sigma_G[, j] <- sigma_G[, j]**abs(a)
  }

  # Draw positively correlated features.
  G <- rmvnorm(n=n, mean=rep(0, n_G_attr), sigma=sigma_G)
  #G <- replicate(n=n_G_attr, rnorm(n=n, mean=0, sd=.15))

  # Standard-normal distribution
  G <- rmvnorm(n=n, mean=rep(0, n_G_attr), sigma=diag(rep(1, n_G_attr)))
  # Data-generating process

  # Degree of sparseness in attributes explaining treatment only
  # Control model complexity of D. With many F, less variance of D
  # is explained by G. This is an important aspect of the DGP.
  sparsity_rate_F <- .9
  beta_F <- rep(beta_F_size, dim(F)[2])
  zero_coefs_F <- sample(1:dim(F)[2], size=sparsity_rate_F*dim(F)[2])
  beta_F[zero_coefs_F] <- 0


  # "Unconfoundedness rate" :
  # Sparsity in both models explaining treatment and outcome,
  # i.e. the same geographical features are important to predict D and Y.
  # However, different effect sizes are allowed.
  sparsity_rate_GDY <- unconfoundedness_rate
  zero_coefs_GDY <- sample(1:dim(G)[2], size=sparsity_rate_GDY*dim(G)[2])

  beta_GD <- rep(beta_GD_size, dim(G)[2])
  beta_GD[zero_coefs_GDY] <- 0

  beta_GY <- rep(beta_GY_size, dim(G)[2])
  beta_GY[zero_coefs_GDY] <- 0


  # Degree of sparseness in housing attributes explaining house prices
  sparsity_rate_H <- .9
  beta_H <- rep(1, dim(H)[2])
  zero_coefs_H <- sample(1:dim(H)[2], size=sparsity_rate_H*dim(H)[2])
  beta_H[zero_coefs_H] <- 0

  #TODO TODO Get this code block organized!
  # Collect information which features' coefficients are set to zero
  sparsity_identifier_F <- rep(0, dim(F)[2])
  sparsity_identifier_G <- rep(0, dim(G)[2])
  sparsity_identifier_H <- rep(0, dim(H)[2])
  sparsity_identifier_F[zero_coefs_F] <- 1
  sparsity_identifier_G[zero_coefs_GDY] <- 1
  sparsity_identifier_H[zero_coefs_H] <- 1
  sparsity_identifier <- append(sparsity_identifier_F, sparsity_identifier_G)
  sparsity_identifier <- append(sparsity_identifier, sparsity_identifier_H)
  sparsity_identifier <- append(sparsity_identifier_G, sparsity_identifier_H) #TODO _F
  sparsity_identifier <- sparsity_identifier_G
  # Collect identifier for variables which are explaining either D or Y
  true_covariate_identifier <- seq(from=1, to=(n_F_attr+n_G_attr+n_H_attr), by=1)
  true_covariate_identifier[sparsity_identifier==1] <- NaN
  true_covariate_identifier <- seq(from=1, to=(n_G_attr+n_H_attr), by=1)
  true_covariate_identifier[sparsity_identifier==1] <- NaN
  true_covariate_identifier <- seq(from=1, to=n_G_attr, by=1)
  true_covariate_identifier[sparsity_identifier==1] <- NaN

  # Check whether true model is high-dimensional
  n_assoc_attr_X <- n_F_attr * (1 - sparsity_rate_F) + n_G_attr * (1 - sparsity_rate_GDY) + n_H_attr * (1 - sparsity_rate_H)
  n_assoc_attr_Z <- n_F_attr * (1 - sparsity_rate_F) + n_G_attr * (1 - sparsity_rate_GDY)
  #if (n/2 < n_assoc_attr_X) {
  #  print("Warning: true model is high-dimensional")
  #}
  #if (n/2 < n_assoc_attr_Z) {
  #  print("Warning: true model is high-dimensional")
  #}

  # Treatment effect
  beta_DY <- treatment_effect

  # Add white noise to both models for treatment and outcome
  eps_D <- rnorm(n=n, mean=0, sd=1)
  eps_Y <- rnorm(n=n, mean=0, sd=1)

  noise_D1 <- rbinom(n=n, size=3, prob=.5)
  noise_D2 <- rpois(n=n, lambda=5)

  # Treatment
  # --- linear model
  D <- 0 + G %*% beta_GD + eps_D
  #D <- 0 + G %*% beta_GD + eps_D + noise_D1 + noise_D2
  #D <- 0 + G %*% beta_GD + F %*% beta_F + eps_D
  # --- logistic model
  #proba_binom_model <- exp(0 + G %*% beta_GD + eps_D) / (1 + exp(0 + G %*% beta_GD + eps_D))
  #D <- cbind(rbinom(n=n, size=1, prob=proba_binom_model))

  # Outcome
  #y <- 0 + D %*% beta_DY + G %*% beta_GY + H %*% beta_H + eps_Y
  y <- 0 + D %*% beta_DY + G %*% beta_GY + eps_Y

  # Add explanatory power outside of model estimated later so that
  # not all variance can be explained.
  noise_y1 <- rbinom(n=n, size=3, prob=.5)
  noise_y2 <- rpois(n=n, lambda=5)
  #y <- 0 + D %*% beta_DY + G %*% beta_GY + H %*% beta_H + eps_Y + noise_y1 + noise_y2
  #TODO



  #TODO
  # Implementation from Belloni et al. (2014)
  n <- n  # 100
  n_G_attr <- n_G_attr  # 200

  beta_DY <- beta_DY  # .5

  # Disturbances
  eps_D <- rnorm(n=n, mean=0, sd=1)
  eps_Y <- rnorm(n=n, mean=0, sd=1)

  # Approximate sparsity
  beta_GD <- 1 / seq(from=1, to=n_G_attr)**2
  beta_GY <- 1 / seq(from=1, to=n_G_attr)**2

  # Regressors
  sigma_G <- matrix(.5, nrow=n_G_attr, ncol=n_G_attr)

  for (j in 1:(n_G_attr)) {
    power <- 1 - j
    a <- seq(from=power, to=(n_G_attr-j), by=1)
    sigma_G[, j] <- sigma_G[, j]**abs(a)
  }

  # Draw positively correlated features.
  #G <- rmvnorm(n=n, mean=rep(0, n_G_attr), sigma=sigma_G)

  # Treatment
  #D <- 0 + G %*% beta_GD + eps_D

  # Outcome
  #y <- 0 + D %*% beta_DY + G %*% beta_GD + eps_Y




  #TODO
  # Implementation of Omitted Variable Bias with Variable Selection
  n <- n  # 500
  n_G_attr <- n_G_attr  # 200

  beta_DY <- beta_DY  # .5

  # Disturbances
  #eps_D <- rnorm(n=n, mean=0, sd=1)
  #eps_Y <- rnorm(n=n, mean=0, sd=1)

  #beta_GD <- rep(1, n_G_attr)
  #beta_GY <- rep(1, n_G_attr)
  #beta_GD[6:n_G_attr] <- 0
  #beta_GY[6:n_G_attr] <- 0

  # Regressors
  #sigma_G <- diag(rep(.05, n_G_attr))

  # Draw features.
  #G <- rmvnorm(n=n, mean=rep(0, n_G_attr), sigma=sigma_G)

  # Treatment
  #D <- 0 + G %*% beta_GD + eps_D

  # Outcome
  #y <- 0 + D %*% beta_DY + G %*% beta_GD + eps_Y
  # Introduce sparseness
  # Set small coefficients to zero
  #a <- round(2/3 * n_G_attr)
  #beta_GD[a:n_G_attr] <- 0
  #beta_GY[a:n_G_attr] <- 0
  #beta_GD[zero_coefs_GDY] <- 0
  #beta_GY[zero_coefs_GDY] <- 0


  # Low dimensions
  #n <- n # 100
  #n_G_attr <- n_G_attr # 1, 2

  #sigma_G <- cbind(c(1, .8), c(.8, 1))
  #X <- rmvnorm(n=n, mean=rep(0, n_G_attr), sigma=sigma_G)
  #G <- c(X[, 1])
  #D <- X[, 2]

  #beta_GD <- .8
  #beta_GY <- .2
  #beta_DY <- 0
  #D <- 0 + beta_GD * G + eps_D
  #y <- 0 + beta_DY * D + beta_GY * G + eps_Y





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
