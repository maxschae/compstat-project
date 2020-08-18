### CompStat project --- Max Schäfer
# Define the data-generating processes


generate_data <- function(dgp=0, n=200, n_F_attr=30, n_G_attr=30, n_H_attr=30,
                          corr_G=0, treatment_effect=.5,
                          beta_GD_size=.5, beta_GY_size=.25, beta_H_size=1, beta_F_size=1,
                          nonzero_controls=NaN) {


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



  # Data-generating process

  # Degree of sparseness in attributes explaining treatment only
  # Control model complexity of D. With many F, less variance of D
  # is explained by G. This is an important aspect of the DGP.
  sparsity_rate_F <- .9
  beta_F <- rep(beta_F_size, dim(F)[2])
  zero_coefs_F <- sample(1:dim(F)[2], size=sparsity_rate_F*dim(F)[2])
  beta_F[zero_coefs_F] <- 0


  # Sparsity in both models explaining treatment and outcome,
  # i.e. the same geographical features are important to predict D and Y.
  # However, different effect sizes are allowed.
  zero_coefs_GDY <- sample(1:n_G_attr, size=(n_G_attr-nonzero_controls))

  beta_GD <- rep(beta_GD_size, n_G_attr)
  beta_GD[zero_coefs_GDY] <- 0

  beta_GY <- rep(beta_GY_size, n_G_attr)
  beta_GY[zero_coefs_GDY] <- 0


  # Degree of sparseness in housing attributes explaining house prices
  sparsity_rate_H <- .9
  beta_H <- rep(1, dim(H)[2])
  zero_coefs_H <- sample(1:dim(H)[2], size=sparsity_rate_H*dim(H)[2])
  beta_H[zero_coefs_H] <- 0



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
  #y <- 0 + D %*% treatment_effect + G %*% beta_GY + H %*% beta_H + eps_Y
  y <- 0 + D %*% treatment_effect + G %*% beta_GY + eps_Y

  # Add explanatory power outside of model estimated later so that
  # not all variance can be explained.
  noise_y1 <- rbinom(n=n, size=3, prob=.5)
  noise_y2 <- rpois(n=n, lambda=5)
  #y <- 0 + D %*% treatment_effect + G %*% beta_GY + H %*% beta_H + eps_Y + noise_y1 + noise_y2
  #TODO



  #TODO
  # Implementation from Belloni et al. (2014)
  # n <- 100
  # n_G_attr <- 200
  # treatment_effect <- .5

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
  #y <- 0 + D %*% treatment_effect + G %*% beta_GD + eps_Y
  # Introduce sparseness
  # Set small coefficients to zero
  #a <- round(2/3 * n_G_attr)
  #beta_GD[a:n_G_attr] <- 0
  #beta_GY[a:n_G_attr] <- 0
  #beta_GD[zero_coefs_GDY] <- 0
  #beta_GY[zero_coefs_GDY] <- 0




  if (dgp == 0) {

    # Specify variance-covariance matrix
    correlation_term <- corr_G
    sigma_G <- matrix(correlation_term, nrow=n_G_attr, ncol=n_G_attr)

    for (j in 1:(n_G_attr)) {
      power <- 1 - j
      a <- seq(from=power, to=(n_G_attr-j), by=1)
      sigma_G[, j] <- sigma_G[, j]**abs(a)
    }

    diag(sigma_G) <- rep(1, n_G_attr)

    # Draw potential controls from multivariate normal distribution.
    G <- rmvnorm(n=n, mean=rep(0, n_G_attr), sigma=sigma_G)

    # Sparsity
    zero_coefs_GDY <- sample(1:n_G_attr, size=(n_G_attr-nonzero_controls))

    beta_GD <- rep(beta_GD_size, n_G_attr)
    beta_GY <- rep(beta_GY_size, n_G_attr)
    beta_GD[zero_coefs_GDY] <- 0
    beta_GY[zero_coefs_GDY] <- 0

    # Disturbances
    eps_D <- rnorm(n=n, mean=0, sd=1)
    eps_Y <- rnorm(n=n, mean=0, sd=1)

    # Treatment
    D <- 0 + G %*% beta_GD + eps_D
    # Outcome
    y <- 0 + D %*% treatment_effect + G %*% beta_GY + eps_Y

    sparsity_identifier <- rep(0, n_G_attr)
    sparsity_identifier[zero_coefs_GDY] <- 1
    true_covariate_identifier <- seq(from=1, to=n_G_attr, by=1)
    true_covariate_identifier[sparsity_identifier==1] <- NaN
  }


  if (dgp == 1) {

    # Specify variance-covariance matrix
    correlation_term <- corr_G
    sigma_G <- matrix(correlation_term, nrow=n_G_attr, ncol=n_G_attr)

    for (j in 1:(n_G_attr)) {
      power <- 1 - j
      a <- seq(from=power, to=(n_G_attr-j), by=1)
      sigma_G[, j] <- sigma_G[, j]**abs(a)
    }

    sigma_F <- diag(rep(1, n_F_attr))
    diag(sigma_G) <- rep(1, n_G_attr)
    sigma_H <- diag(rep(1, n_H_attr))

    # Draw potential controls from multivariate normal distribution.
    F <- rmvnorm(n=n, mean=rep(0, n_F_attr), sigma=sigma_F)
    G <- rmvnorm(n=n, mean=rep(0, n_G_attr), sigma=sigma_G)
    H <- rmvnorm(n=n, mean=rep(0, n_H_attr), sigma=sigma_H)

    # Sparsity
    zero_coefs_F <- sample(1:n_F_attr, size=(n_F_attr-nonzero_controls))
    zero_coefs_GDY <- sample(1:n_G_attr, size=(n_G_attr-nonzero_controls))
    zero_coefs_H <- sample(1:n_H_attr, size=(n_H_attr-nonzero_controls))

    # Coefficients
    beta_F <- rep(beta_F_size, n_F_attr)
    beta_GD <- rep(beta_GD_size, n_G_attr)
    beta_GY <- rep(beta_GY_size, n_G_attr)
    beta_H <- rep(beta_H_size, n_H_attr)

    beta_F[zero_coefs_F] <- 0
    beta_GD[zero_coefs_GDY] <- 0
    beta_GY[zero_coefs_GDY] <- 0
    beta_H[zero_coefs_H] <- 0

    # Disturbances
    eps_D <- rnorm(n=n, mean=0, sd=1)
    eps_Y <- rnorm(n=n, mean=0, sd=1)

    # Treatment
    D <- 0 + F %*% beta_F + G %*% beta_GD + eps_D
    # Outcome
    y <- 0 + D %*% treatment_effect + G %*% beta_GY + H %*% beta_H + eps_Y



    sparsity_identifier_F <- rep(0, n_F_attr)
    sparsity_identifier_G <- rep(0, n_G_attr)
    sparsity_identifier_H <- rep(0, n_H_attr)
    sparsity_identifier_F[zero_coefs_F] <- 1
    sparsity_identifier_G[zero_coefs_GDY] <- 1
    sparsity_identifier_H[zero_coefs_H] <- 1
    sparsity_identifier <- append(sparsity_identifier_F,
                                  sparsity_identifier_G)
    sparsity_identifier <- append(sparsity_identifier, sparsity_identifier_H)

    true_covariate_identifier <- seq(from=1, to=(n_F_attr+n_G_attr+n_H_attr), by=1)
    true_covariate_identifier[sparsity_identifier==1] <- NaN
  }





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
