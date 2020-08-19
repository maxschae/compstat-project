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



generate_data_houseprices <- function(n=200, n_F_attr=30, n_G_attr=100, n_H_attr=100, corr_G=0) {


  # Data-generating process inspired by [PAPER]



  # Potential covariates
  F_1 <- rchisq(n=n, df=1) * 2500
  beta_F <- c(.25)

  # Potential confounding variables
  # 1x1-km cluster : n / k
  G_1 <- rchisq(n=n, df=3) * 2000 + 15000   # Purchasing power
  G_2 <- rchisq(n=n, df=2) * 1000   # Inhabitants
  G_3 <- abs(rchisq(n=n, df=6) * 80 - 40)  # Number of buildings
  G_4 <- rbeta(n=n, 3, 80) * 1000  # Distance to city center (in km)
  G <- cbind(G_1, G_2, G_3, G_4)

  # Coefficients
  beta_GD <- c(1, -.5, .1, -1)
  beta_GY <- c(.0382, .031, -.00002, -.0042)

  # House characteristics
  H_1 <- runif(n=n, min=1945, max=2016)   # Year of construction
  H_2 <- rchisq(n=n, df=3) * 25 + 80  # Living space
  H_3 <- rchisq(n=n, df=5) * 140  # Lot size
  H_4 <- rbinom(n=n, size=12, prob=.3) + 2  # Number of rooms

  H <- cbind(H_1, H_2, H_3, H_4)

  # House-type indicator variables
  H_5d1 <- rep(0, n)
  H_5d2 <- rep(0, n)
  #H_5d3 <- rep(0, n)   # Leave out one house-type to avoid perfect multicollinearity
  H_5d4 <- rep(0, n)
  H_5d5 <- rep(0, n)
  H_5d6 <- rep(0, n)
  H_5d7 <- rep(0, n)
  H_5d8 <- rep(0, n)

  H_5d1[seq(from=0, to=.55*n, by=1)] <- 1
  H_5d2[seq(from=round(.55*n), to=round((.55+.17)*n), by=1)] <- 1
  #H_5d3[seq(from=round((.55+.17)*n), to=round((.55+.17+.08)*n), by=1)] <- 1
  H_5d4[seq(from=round((.55+.17+.08)*n), to=round((.55+.17+.08+.04)*n), by=1)] <- 1
  H_5d5[seq(from=round((.55+.17+.08+.04)*n), to=round((.55+.17+.08+.04+.06)*n), by=1)] <- 1
  H_5d6[seq(from=round((.55+.17+.08+.04+.06)*n), to=round((.55+.17+.08+.04+.06+.04)*n), by=1)] <- 1
  H_5d7[seq(from=round((.55+.17+.08+.04+.06+.04)*n), to=round((.55+.17+.08+.04+.06+.04)*n), by=1)] <- 1
  H_5d8[seq(from=round((.55+.17+.08+.04+.06+.04+.03)*n), to=n, by=1)] <- 1

  H_5d <- cbind(H_5d1, H_5d2, H_5d4, H_5d5, H_5d6, H_5d7, H_5d8)

  # Coefficients
  beta_H <- c(.00453, .0041, .0122, -.0104)
  beta_H5d <- c(.026, -.0601, -.153, -.149, -.0838, .0383, .214)


  # Add some noise two both stages
  eps_D <- rnorm(n=n, mean=0, sd=100)
  eps_Y <- rnorm(n=n, mean=0, sd=1000)

  # Treatment (distance to next wind turbine in km)
  # DGP-0
  D <- (-12000 + G %*% beta_GD + eps_D) #/ 1000

  # DGP-1
  #D <- (0 + beta_G %*% cbind(G_1, G_2, G_3, G_4) + beta_F %*% F_1 + eps_D) / 1000

  beta_DY <- -.004

  # Outcome
  # In the paper, the house price equation is expressed in log house prices
  # and 'regular' units on the right-hand-side. Here, house prices are in EUR
  # and coefficients are scaled down up by roughly one-percent of the median house price.
  c <- .01*175000
  treatment_effect <- beta_DY*c
  y <- -1200000 + D %*% treatment_effect + G %*% beta_GY*c + H %*% beta_H*c + H_5d %*% beta_H5d*c + eps_Y



  # Construct geographical features and house attributes
  # which do not explain anything
  #zero_coefs <- n_F_attr + n_G_attr + n_H_attr - dim(H)[2] - dim(H_5d)[2] - dim(G)[2] - dim(F)[2]
  n_G_attr_zero <- n_G_attr - dim(G)[2]
  n_H_attr_zero <- n_H_attr - dim(H)[2] - dim(H_5d)[2]


  # Specify variance-covariance matrix
  correlation_term <- .5
  sigma_G_zero <- matrix(correlation_term, nrow=n_G_attr_zero, ncol=n_G_attr_zero)

  # Covariances must be smaller than variances for covariance matrix to be positive semi-definite
  for (j in 1:(n_G_attr_zero)) {
    power <- 1 - j
    a <- seq(from=power, to=(n_G_attr_zero-j), by=1)
    sigma_G_zero[, j] <- 100*sigma_G_zero[, j]**abs(a)
  }

  #sigma_F <- diag(rep(1, n_F_attr_zero))
  diag(sigma_G_zero) <- rep(100, n_G_attr_zero)
  sigma_H_zero <- diag(rep(100, n_H_attr_zero))

  # Draw potential controls from multivariate normal distribution.
  #F <- rmvnorm(n=n, mean=rep(0, n_F_attr), sigma=sigma_F_zero)
  G_zero <- rmvnorm(n=n, mean=rep(0, n_G_attr_zero), sigma=sigma_G_zero)
  H_zero <- rmvnorm(n=n, mean=rep(0, n_H_attr_zero), sigma=sigma_H_zero)


  # Alternatives
  #G_1 <- rchisq(n=n, df=5) * 4200 + 5000
  #G_1 <- rchisq(n=n, df=2) * 4200 + 13000
  #G_1 <- rchisq(n=n, df=10) * 1100 + 15000
  #G_2 <- abs(rbeta(n=n, 3, 50) * 50000 - 1000) + 1
  #G_4 <- abs(rbeta(n=n, 2, 15) * 260 - 5)



  #sparsity_identifier_F <- rep(0, n_F_attr)
  sparsity_identifier_G <- rep(0, n_G_attr_zero)
  sparsity_identifier_H <- rep(0, n_H_attr_zero)
  #sparsity_identifier_F[zero_coefs_F] <- 1
  sparsity_identifier_G[(dim(G)[2]+1):n_G_attr] <- 1
  sparsity_identifier_H[(dim(H)[2]+dim(H_5d)[2]+1):n_H_attr] <- 1
  #sparsity_identifier <- append(sparsity_identifier_F,
  #                              sparsity_identifier_G)
  sparsity_identifier <- append(sparsity_identifier_G, sparsity_identifier_H)

  true_covariate_identifier <- seq(from=1, to=(n_G_attr+n_H_attr), by=1)
  true_covariate_identifier[sparsity_identifier==1] <- NaN

  # Create vector with covariate names.
  #colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
  colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
  colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))

  # Required for first selection step
  #colnames_dataset <- c("y", "D", colname_F, colname_G, colname_H)
  colnames_dataset <- c("y", "D", colname_G, colname_H)

  #data <- data.frame(cbind(y, D, F, G, G_zero, H, H_5d, H_zero))
  data <- data.frame(cbind(y, D, G, G_zero, H, H_5d, H_zero))
  names(data) <- colnames_dataset
  #return(list(y=y, D=D, F=F, G=G, H=H, true_covariate_identifier=true_covariate_identifier))
  return(list(data=data, true_covariate_identifier=true_covariate_identifier))
}


generate_data_houseprices_dgp0 <- function(n=200, n_G_attr=100, corr_G=0, treatment_effect=-7, beta_GY_inflator=1) {


  # Data-generating process inspired by [PAPER]


  # Potential confounding variables
  # 1x1-km cluster : n / k
  G_1 <- rchisq(n=n, df=3) * 2000 + 15000   # Purchasing power
  G_2 <- rchisq(n=n, df=2) * 1000   # Inhabitants
  G_3 <- abs(rchisq(n=n, df=6) * 80 - 40)  # Number of buildings
  G_4 <- rbeta(n=n, 3, 80) * 1000  # Distance to city center (in km)
  G <- cbind(G_1, G_2, G_3, G_4)

  # Coefficients
  beta_GD <- c(1, -.5, .1, -1)
  beta_GY <- c(.0382, .031, -.00002, -.0042)

  beta_GY <- beta_GY * beta_GY_inflator

  # Add some noise two both stages
  eps_D <- rnorm(n=n, mean=0, sd=100)
  eps_Y <- rnorm(n=n, mean=0, sd=1000)

  # Treatment (distance to next wind turbine in km)
  # DGP-0
  D <- (-12000 + G %*% beta_GD + eps_D) #/ 1000

  # DGP-1
  #D <- (0 + beta_G %*% cbind(G_1, G_2, G_3, G_4) + beta_F %*% F_1 + eps_D) / 1000



  # Outcome
  # In the paper, the house price equation is expressed in log house prices
  # and 'regular' units on the right-hand-side. Here, house prices are in EUR
  # and coefficients are scaled down up by roughly one-percent of the median house price.
  c <- .01*175000
  beta_DY <- -.004
  #treatment_effect <- beta_DY*c
  treatment_effect <- treatment_effect
  y <- -1200000 + D %*% treatment_effect + G %*% beta_GY*c + eps_Y



  # Construct geographical features and house attributes
  # which do not explain anything
  n_G_attr_zero <- n_G_attr - dim(G)[2]


  # Specify variance-covariance matrix
  correlation_term <- .5
  sigma_G_zero <- matrix(correlation_term, nrow=n_G_attr_zero, ncol=n_G_attr_zero)

  # Covariances must be smaller than variances for covariance matrix to be positive semi-definite
  for (j in 1:(n_G_attr_zero)) {
    power <- 1 - j
    a <- seq(from=power, to=(n_G_attr_zero-j), by=1)
    sigma_G_zero[, j] <- 100*sigma_G_zero[, j]**abs(a)
  }

  #sigma_F <- diag(rep(1, n_F_attr_zero))
  diag(sigma_G_zero) <- rep(100, n_G_attr_zero)

  # Draw potential controls from multivariate normal distribution.
  #F <- rmvnorm(n=n, mean=rep(0, n_F_attr), sigma=sigma_F_zero)
  G_zero <- rmvnorm(n=n, mean=rep(0, n_G_attr_zero), sigma=sigma_G_zero)


  sparsity_identifier_G <- rep(0, n_G_attr_zero)
  sparsity_identifier_G[(dim(G)[2]+1):n_G_attr] <- 1

  sparsity_identifier <- sparsity_identifier_G

  true_covariate_identifier <- seq(from=1, to=n_G_attr, by=1)
  true_covariate_identifier[sparsity_identifier==1] <- NaN

  # Create vector with covariate names.
  colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))

  # Required for first selection step
  colnames_dataset <- c("y", "D", colname_G)

  #data <- data.frame(cbind(y, D, F, G, G_zero, H, H_5d, H_zero))
  data <- data.frame(cbind(y, D, G, G_zero))
  names(data) <- colnames_dataset
  #return(list(y=y, D=D, F=F, G=G, H=H, true_covariate_identifier=true_covariate_identifier))
  return(list(data=data, true_covariate_identifier=true_covariate_identifier))
}












#
