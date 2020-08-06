### Data-generating process
rm(list=ls())

set.seed(1234)

# install.packages("glmnet")
library(glmnet)
library(mvtnorm)

library(stringr)

library(ggplot2)
#install.packages("randomcoloR")
library(randomcoloR)

# Own modules
#source("estimation.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/estimation.R")



generate_data <- function(n=200, n_G_attr, n_H_attr) {


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


  # Treatment effect
  beta_DY <- -1


  # Treatment
  # --- linear model
  D <- 0 + G %*% beta_GD + eps_D
  # --- logistic model
  #proba_binom_model <- exp(0 + G %*% beta_GD + eps_D) / (1 + exp(0 + G %*% beta_GD + eps_D))
  #D <- cbind(rbinom(n=n, size=1, prob=proba_binom_model))

  # Outcome
  y <- 0 + D %*% beta_DY + G %*% beta_GY + H %*% beta_H + eps_Y


  return(list(y=y, D=D, G=G, H=H))
}


### Main ---

n <- 200
n_G_attr <- 10
n_H_attr <- 15


# Create vector with covariate names.
colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))
colnames_DG <- c("D", colname_G)
colnames_GH <- c(colname_G, colname_H)
colnames_DGH <- c("D", colname_G, colname_H)



data <- generate_data(n=n, n_G_attr=n_G_attr, n_H_attr=n_H_attr)

y <- data$y
D <- data$D
G <- data$G
H <- data$H


# Regression results
lm(y ~ D + G + H)
lm(y ~ D + H)
lm(y ~ D)

lm(D ~ G + H)
lm(D ~ G)


beta_hat <- least_squares_estimator(y=y, X=cbind(D, G, H))
beta_hat <- least_squares_estimator(y=y, X=cbind(D, H))
beta_hat <- least_squares_estimator(y=y, X=D)



### LASSO


# Data
data_train <- generate_data(n=n, n_G_attr=n_G_attr, n_H_attr=n_H_attr)

y_train <- data_train$y
D_train <- data_train$D
G_train <- data_train$G
H_train <- data_train$H

X_train <- cbind(D_train, G_train, H_train)


# Double-selection method
lambda_grid <- seq(0, 2, by=.025)
# 1. Use Lasso to select covariates based their association with outcome Y
lasso_one <- glmnet(X_train, y_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

# Extract estimates from lasso model
beta_hats_lasso_one <- t(matrix(coef(lasso_one, s=lambda_grid)[-1,], nrow=dim(X_train)[2], ncol=length(lambda_grid)))



#TODO Quick inspection
#plot(lasso_one)
#plot(lasso_one, label=TRUE)
#print(lasso_one)
#coef(lasso_one, s=.1)


# Using METRIC-optimal lambda
lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE)

selected_covars_one <- coef(lasso_one_cv, s="lambda.min")[-1]
#lambda_cv <- lasso_one_cv$lambda_min
#selection_one <- coef(lasso_one, s=lambda_cv)


# Create vector that indicates selected covariates by lasso_one.
column_indicator_covars_one <- seq(from=1, to=length(selected_covars_one), by=1)
column_indicator_covars_one <- t(cbind(column_indicator_covars_one))
# Name covariates.
names(column_indicator_covars_one) <- colnames_DGH
# Keep indicators of selected covariates.
column_indicator_covars_one[selected_covars_one==0] <- NaN
column_indicator_covars_one <- column_indicator_covars_one[!is.na(column_indicator_covars_one)]

# Drop covariates from data that were excluded by Lasso.
X_train_select_one <- X_train[1:n, column_indicator_covars_one]





# 2. Use Lasso to select covariates given their association with treatment D
#TODO cannot use logistic model then?? or use extension of lasso with logit...
# type.measure class gives misclassification error

#TODO use X_train or G_train, since H_train should intuitively have no impact
# for wind turbine --- but perhaps for other treatments as school placements?
#lasso_two <- glmnet(G_train, D_train, alpha=1, family="binomial", intercept=FALSE)
#lasso_two_cv <- cv.glmnet(G_train, D_train, alpha=1, family="binomial", type.measure="class", intercept=FALSE)
GH_train <- cbind(G_train, H_train)

lasso_two <- glmnet(G_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
lasso_two <- glmnet(GH_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

# Extract estimates from lasso model
#beta_hats_lasso_two <- t(matrix(coef(lasso_two, s=lambda_grid)[-1,], nrow=dim(G_train)[2], ncol=length(lambda_grid)))
beta_hats_lasso_two <- t(matrix(coef(lasso_two, s=lambda_grid)[-1,], nrow=dim(GH_train)[2], ncol=length(lambda_grid)))








# Plot all lasso coefficients with increasing lambda
lasso_coef_shrink_plot <- function(df) {

  cols <- names(df)[1:(length(names(df))-1)]

  # Initialize plot.
  plot <- ggplot(data=df, aes(x=lambda_grid))
  # Loop through each beta and plot its value with increasing penalty.
  for (name in cols) {
    random_transparency <- runif(n=1, min=.5, max=1) #
    if (startsWith(name, "D")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#43f7d0")
    }
    else if (startsWith(name, "G")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#f74376", alpha=random_transparency)
    }
    else if (startsWith(name, "H")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#6143f7", alpha=random_transparency)
    }
  }
  plot <- plot +
          geom_line(aes(y=0), lty="dashed") +
          ylab("Lasso coefficients") +
          xlab("Lambda")

  plot
}

# 1.
df <- data.frame(cbind(beta_hats_lasso_one, lambda_grid))
names(df) <- colnames_DGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

lasso_coef_shrink_plot(df=df)


# 2.
df <- data.frame(cbind(beta_hats_lasso_two, lambda_grid))
names(df) <- colnames_GH
#names(df) <- colnames_DG
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

lasso_coef_shrink_plot(df=df)






# OLS with union of X_select_one and X_select_two


























### Code snippets

# Ensure correlations are positive
#A <- matrix(rnorm(n=H_no_attr**2, mean=0, sd=2), ncol=H_no_attr)
#sigma <- t(A) %*% A
#diag(sigma) <- diag(sigma) * 10
#sigma <- abs(sigma)

#x <- rbinom(n=10, size=100, prob=.15)
#print(x)


# Initialize plot.
#plot <- ggplot(data=df, aes(x=lambda_grid))
# Loop through each beta and plot its value with increasing penalty.
#for (name in cols) {
#  plot <- plot +
#          geom_line(aes_string(y=name), size=1, col=randomColor(count=1, luminosity="light"))
#}
#plot <- plot +
#        geom_line(aes(y=0), lty="dashed") +
#        ylab("Lasso coefficients") +
#        xlab("Lambda")
#plot



#
