### CompStat project --- Max Schäfer
# Define the data-generating processes


generate_data_A <- function(n=400, n_G_attr=100, corr_G=0, treatment_effect=.25,
                            beta_GD_size=.25, beta_GY_size=.1,
                            nonzero_controls=NaN) {

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

  # Direct effect sizes
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

  # Create vector with covariate names.
  colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))

  # Required for first selection step
  colnames_dataset <- c("y", "D", colname_G)
  data <- data.frame(cbind(y, D, G))
  names(data) <- colnames_dataset

  return(list(data=data, true_covariate_identifier=true_covariate_identifier))
}




generate_data_houseprices <- function(n=200, n_F_attr=100, n_G_attr=100, n_H_attr=100,
                                      beta_GY_inflator=1, beta_GD_inflator=1, corr_G=.1) {

  # Data-generating process inspired by Frondel et al. (2019)


  # Potential covariates
  F_1 <- rchisq(n=n, df=1) * 2500

  F_random <- rnorm(n=n, mean=0, sd=1)

  F <- cbind(F_1, F_random)
  beta_F <- c(0,  0)

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

  beta_GD <- beta_GD * beta_GD_inflator
  # Treatment (distance to next wind turbine in km)
  # DGP-1
  D <- (-12000 + G %*% beta_GD + F %*% beta_F + eps_D) #/ 1000


  beta_DY <- -.004

  # Outcome
  # In the paper, the house price equation is expressed in log house prices
  # and 'regular' units on the right-hand-side. Here, house prices are in EUR
  # and coefficients are scaled up by roughly one-percent of the median house price.
  c <- .01*175000
  treatment_effect <- beta_DY*c
  beta_GY <- beta_GY * beta_GY_inflator

  y <- -1200000 + D %*% treatment_effect + G %*% beta_GY*c + H %*% beta_H*c + H_5d %*% beta_H5d*c + eps_Y



  # Construct geographical features and house attributes
  # which do not explain anything
  n_F_attr_zero <- n_F_attr - dim(F)[2]
  n_G_attr_zero <- n_G_attr - dim(G)[2]
  n_H_attr_zero <- n_H_attr - dim(H)[2] - dim(H_5d)[2]

  # Specify variance-covariance matrix
  correlation_term <- 0
  sigma_G_zero <- matrix(correlation_term, nrow=n_G_attr_zero, ncol=n_G_attr_zero)

  # Covariances must be smaller than variances for covariance matrix to be positive semi-definite
  for (j in 1:(n_G_attr_zero)) {
    power <- 1 - j
    a <- seq(from=power, to=(n_G_attr_zero-j), by=1)
    sigma_G_zero[, j] <- 100*sigma_G_zero[, j]**abs(a)
  }

  sigma_F_zero <- diag(rep(100, n_F_attr_zero))
  diag(sigma_G_zero) <- rep(100, n_G_attr_zero)
  sigma_H_zero <- diag(rep(100, n_H_attr_zero))

  # Draw potential controls from multivariate normal distribution.
  F_zero <- rmvnorm(n=n, mean=rep(0, n_F_attr_zero), sigma=sigma_F_zero)
  G_zero <- rmvnorm(n=n, mean=rep(0, n_G_attr_zero), sigma=sigma_G_zero)
  H_zero <- rmvnorm(n=n, mean=rep(0, n_H_attr_zero), sigma=sigma_H_zero)


  sparsity_identifier_F <- rep(0, n_F_attr_zero)
  sparsity_identifier_G <- rep(0, n_G_attr_zero)
  sparsity_identifier_H <- rep(0, n_H_attr_zero)
  sparsity_identifier_F[(dim(G)[2]+1):n_F_attr] <- 1
  sparsity_identifier_G[(dim(G)[2]+1):n_G_attr] <- 1
  sparsity_identifier_H[(dim(H)[2]+dim(H_5d)[2]+1):n_H_attr] <- 1
  sparsity_identifier <- append(sparsity_identifier_F, sparsity_identifier_G)
  sparsity_identifier <- append(sparsity_identifier, sparsity_identifier_H)

  true_covariate_identifier <- seq(from=1, to=(n_F_attr+n_G_attr+n_H_attr), by=1)
  true_covariate_identifier[sparsity_identifier==1] <- NaN

  # Create vector with covariate names.
  colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
  colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
  colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))

  # Required for first selection step
  colnames_dataset <- c("y", "D", colname_F, colname_G, colname_H)

  data <- data.frame(cbind(y, D, F, F_zero, G, G_zero, H, H_5d, H_zero))
  names(data) <- colnames_dataset

  return(list(data=data, true_covariate_identifier=true_covariate_identifier))
}


#
