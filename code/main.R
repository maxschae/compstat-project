### CompStat project --- Max Schäfer

rm(list=ls())


#install.packages("glmnet")
#install.packages("randomcoloR")

library(glmnet)
library(mvtnorm)

library(stringr)

library(ggplot2)
library(randomcoloR)

# Own modules
#source("estimation.R")
#source("data_generating_process.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/estimation.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/data_generating_process.R")



### Main ---
set.seed(1234)


n <- 400
n_F_attr <- 70
n_G_attr <- 90
n_H_attr <- 90
treatment_effect <- -1
unconfoundedness_rate <- .5

# Housekeeping of feature names and identifiers for identification.
# Create vector with covariate names.
colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))

# Required for first selection step
colnames_DGH <- c("D", colname_G, colname_H)
colnames_DFGH <- c("D", colname_F, colname_G, colname_H)
colnames_dataset <- c("y", "D", colname_F, colname_G, colname_H)

# Required for second selection step
colnames_GH <- c(colname_G, colname_H)
colnames_FGH <- c(colname_F, colname_G, colname_H)

colnumbers_GH <- seq(from=1, to=(n_G_attr+n_H_attr), by=1)
colnumbers_FGH <- seq(from=1, to=(n_F_attr+n_G_attr+n_H_attr), by=1)




### Double-selection method

simple_select_covariates <- function(data_train) {


  y_train <- data_train$y
  D_train <- data_train$D
  F_train <- data_train[colname_F]
  G_train <- data_train[colname_G]
  H_train <- data_train[colname_H]


  # Assemble different regressor matrices
  # Case with F_
  X_train <- data.matrix(cbind(D_train, F_train, G_train, H_train))
  Z_train <- data.matrix(cbind(F_train, G_train, H_train))
  # Case without F_
  X_train <- data.matrix(cbind(D_train, G_train, H_train))
  Z_train <- data.matrix(cbind(G_train, H_train))

  # Simple selection method
  # Off-the-shelve MSE-optimal Lasso; outcome regressed on entire feature matrix

  # Use Lasso to select covariates given their association with outcome Y
  lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.min")[-1])
  # Exclude coef. of treatment since we are interested in selecting covariates
  lasso_coefs_one <- lasso_coefs_one[2:dim(X_train)[2], ]

  # Create vector that indicates selected covariates by lasso_one_cv
  selected_covars_one <- seq(from=1, to=length(lasso_coefs_one), by=1)
  selected_covars_one <- rbind(selected_covars_one)
  # Name covariates.
  names(selected_covars_one) <- colnames_GH
  #names(selected_covars_one) <- colnames_FGH
  # Keep indicators of selected covariates.
  selected_covars_one[lasso_coefs_one==0] <- NaN

  return(selected_covars_one)
}



double_select_covariates <- function(data_train) {

  y_train <- data_train$y
  D_train <- data_train$D
  F_train <- data_train[colname_F]
  G_train <- data_train[colname_G]
  H_train <- data_train[colname_H]


  # Assemble different regressor matrices
  # Case with F_
  X_train <- data.matrix(cbind(D_train, F_train, G_train, H_train))
  Z_train <- data.matrix(cbind(F_train, G_train, H_train))
  # Case without F_
  X_train <- data.matrix(cbind(D_train, G_train, H_train))
  Z_train <- data.matrix(cbind(G_train, H_train))

  #GH_train <- data.matrix(cbind(G_train, H_train))
  #FGH_train <- data.matrix(cbind(F_train, G_train, H_train))


  # Double-selection method

  # 1. Use Lasso to select covariates given their association with outcome Y
  selected_covars_one <- simple_select_covariates(data_train=data_train)

  # 2. Use Lasso to select covariates given their association with tratment D
  lasso_two_cv <- cv.glmnet(Z_train, D_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lasso_coefs_two <- cbind(coef(lasso_two_cv, s="lambda.min")[-1])

  # Create vector that indicates selected covariates by lasso_two_cv
  selected_covars_two <- seq(from=1, to=length(lasso_coefs_two), by=1)
  selected_covars_two <- t(cbind(selected_covars_two))
  # Name covariates
  names(selected_covars_two) <- colnames_GH
  #names(selected_covars_two) <- colnames_FGH
  # Keep indicators of selected covariates.
  selected_covars_two[lasso_coefs_two==0] <- NaN
  #selected_covars_two <- selected_covars_two[!is.na(selected_covars_two)]


  # Exclude covariates only if both lasso_one and lasso_two suggest it
  double_selected_identifier <- rbind(selected_covars_one,
                                      selected_covars_two,
                                      colnumbers_GH)
                                      #colnumbers_FGH

  double_selected_identifier[3, colSums(is.na(double_selected_identifier)) == 2] <- NaN
  double_selected_identifier <- double_selected_identifier[3, ]

  return(double_selected_identifier)
}


###

compute_confounder_bias <- function(n, n_F_attr, n_G_attr, n_H_attr, treatment_effect,
                                    unconfoundedness_rate) {


  # Draw data once and split into train and test (*)
  data_set <- generate_data(n=n,
                        n_F_attr=n_F_attr,
                        n_G_attr=n_G_attr,
                        n_H_attr=n_H_attr,
                        treatment_effect=treatment_effect,
                        unconfoundedness_rate=unconfoundedness_rate)

  data <- data_set$data

  # Split data into training and test sets
  data_train <- data[1:(dim(data)[1]/2), ]
  data_test <- data[(dim(data)[1]/2+1):dim(data)[1], ]

  true_covariate_identifier <- data_set$true_covariate_identifier

  # *Selection section*
  # Select covariates according to the double-selection method
  simple_selected_identifier <- simple_select_covariates(data_train=data_train)
  double_selected_identifier <- double_select_covariates(data_train=data_train)


  covariate_identifier <- data.frame(rbind(true_covariate_identifier,
                                           simple_selected_identifier,
                                           double_selected_identifier))

  #names(covariate_identifier) <- colnames_FGH
  names(covariate_identifier) <- colnames_GH

  #covariate_identifier_F <- covariate_identifier[colname_F] #TODO F_
  covariate_identifier_G <- covariate_identifier[colname_G]
  covariate_identifier_H <- covariate_identifier[colname_H]



  #TODO check whether treatment was also dropped?

  # Be informed whether variabels were falsely dropped
  #TODO False positive and false negative
  #misexcl_F_rate <- sum(is.na(covariate_identifier_F[1,] != covariate_identifier_F[2,])) / n_F_attr
  misexcl_G_rate_simple <- sum(is.na(covariate_identifier_G[1,] != covariate_identifier_G[2,])) / n_G_attr
  misexcl_H_rate_simple <- sum(is.na(covariate_identifier_H[1,] != covariate_identifier_H[2,])) / n_H_attr

  misexcl_G_rate_double <- sum(is.na(covariate_identifier_G[1,] != covariate_identifier_G[3,])) / n_G_attr
  misexcl_H_rate_double <- sum(is.na(covariate_identifier_H[1,] != covariate_identifier_H[3,])) / n_H_attr
  #TODO check this

  # The setup here is labeled as follows:
  # . 'positive' refers to state where covariate has a zero effect.
  # . 'negative' refers to state where covariate has a non-zero effect.
  # Thus, true-positive refers to state where covariate has a zero effect
  # and the selection methods excludes covariate.
  # True-negative refers to state where covariate has a non-zero effect
  # and the selection method does not exclude covariate.
  # False-positive : covariate has a zero effect and is not excluded.
  # False-negative : covariate has non-zero effect and is not excluded.

  tp_rate_G <- 1


  # Set in context to sparsity rates
  #misexcl_F_rate <- sum(is.na(covariate_identifier_F[1,] != covariate_identifier_F[2,])) / n_F_attr / sparsity_rate_F
  #misexcl_G_rate <- sum(is.na(covariate_identifier_G[1,] != covariate_identifier_G[2,])) / n_G_attr / sparsity_rate_GDY
  #misexcl_H_rate <- sum(is.na(covariate_identifier_H[1,] != covariate_identifier_H[2,])) / n_H_attr / sparsity_rate_H


  # Compute metrics of interest (see below) with test data.
  y_test <- data_test$y
  D_test <- data_test$D
  F_test <- data_test[colname_F]
  G_test <- data_test[colname_G]
  H_test <- data_test[colname_H]

  # Assemble different regressor matrices
  #X_test <- data.matrix(cbind(D_test, F_test, G_test, H_test))
  X_test <- data.matrix(cbind(D_test, G_test, H_test))
  Z_test <- data.matrix(cbind(G_test, H_test))
  #Z_test <- data.matrix(cbind(F_test, G_test, H_test))

  #Z_test_selected <- Z_test[, double_selected_identifier[!is.na(double_selected_identifier)]]
  Z_test_selected <- Z_test[, simple_selected_identifier[!is.na(simple_selected_identifier)]]

  X_test_selected <- data.matrix(cbind(D_test, Z_test_selected))

  # Estimate treatment effect by standard OLS provided with double-selected covariates
  lm_result <- lm(y_test ~ X_test_selected)
  # Collect metrics of interest
  confounder_bias <- lm_result$coefficients[2] - treatment_effect
  test_MSE <- (1/length(lm_result$residuals)) * sum((lm_result$residuals)**2)

  return(list(confounder_bias=confounder_bias,
              test_MSE=test_MSE,
              misexcl_G_rate_simple=misexcl_G_rate_simple,
              misexcl_H_rate_simple=misexcl_H_rate_simple))
  # (*) Since variables (within a class like G or H) are
  # excluded from the true model at random, drawing *one* dataset which
  # is then split into train and test data ensures that the true association
  # of variables is equal across train and test data.
  # E.g. if train and test data are drawn sequentially, "G_1" might be
  # associated in the train but not in the test data, and dropping it from test
  # is fine whereas dropping it in train is not. But then we would
  # reason that "G_1" was falsely dropped albeit it was not. This
  # would also mess up any metric such as the MSE etc.
}


# Main : run simulation
R <- 20


sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                       n_F_attr=n_F_attr,
                                       n_G_attr=n_G_attr,
                                       n_H_attr=n_H_attr,
                                       treatment_effect=treatment_effect,
                                       unconfoundedness_rate=unconfoundedness_rate))

confounder_bias_vec <- unlist(sim_results_vec[1, 1:R])
test_MSE_vec <- unlist(sim_results_vec[2, 1:R])
misexcl_G_rates <- unlist(sim_results_vec[3, 1:R])
misexcl_H_rates <- unlist(sim_results_vec[4, 1:R])



# Evaluate
cat("Mean test MSE is", mean(test_MSE_vec))
cat("Median confounder bias is", median(confounder_bias_vec))

hist(confounder_bias_vec)




# Simulations
# Vary over
unconfoundedness_rate_vec <- seq(from=0, to=1, by=.1)

# Initialize containers
confounder_bias_mat <- matrix(NA, nrow=length(unconfoundedness_rate_vec), ncol=R)
mean_test_MSE_vec <- rep(NA, R)
mean_misexcl_G_rate_simple <- rep(NA, R)
mean_misexcl_H_rate_simple <- rep(NA, R)

i <- 1

for (unconfoundedness_rate in unconfoundedness_rate_vec) {

  sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                         n_F_attr=n_F_attr,
                                         n_G_attr=n_G_attr,
                                         n_H_attr=n_H_attr,
                                         treatment_effect=treatment_effect,
                                         unconfoundedness_rate=unconfoundedness_rate))

  confounder_bias_mat[i, ] <- unlist(sim_results_vec[1, 1:R])
  mean_test_MSE_vec[i] <- mean(unlist(sim_results_vec[2, 1:R]))
  mean_misexcl_G_rate_simple[i] <- mean(unlist(sim_results_vec[3, 1:R]))
  mean_misexcl_H_rate_simple[i] <- mean(unlist(sim_results_vec[4, 1:R]))


  i <- i+1
}


# Evaluation



#TODO do this with replicate function, R times. then: loop / vectorize over it by
# iterating over treatment_effect or numbers of attributes or unconfoundedness rate / intensity etc.
# this forms the key simulation study ---
# loop or vectorize function w.r.t. what should be varied.








# 2. Use Lasso to select covariates given their association with treatment D
#TODO cannot use logistic model then?? or use extension of lasso with logit...
# type.measure class gives misclassification error for cv
#lasso_two_cv <- cv.glmnet(G_train, D_train, alpha=1, family="binomial", type.measure="class", intercept=FALSE)
















###










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
    else if (startsWith(name, "F")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#a89932", alpha=random_transparency)
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


# Plotting section ---



#TODO Do this for multiple runs... quite tricky since its 3-dimensions
# i.e. lambdas, coefs and simulation runs span three dimensions

data_set <- generate_data(n=n,
                          n_F_attr=n_F_attr,
                          n_G_attr=n_G_attr,
                          n_H_attr=n_H_attr,
                          treatment_effect=treatment_effect,
                          unconfoundedness_rate=unconfoundedness_rate)

data <- data_set$data
data_train <- data[1:dim(data)[1]/2, ]

# Draw training data
y_train <- data_train$y
D_train <- data_train$D
F_train <- data_train[colname_F]
G_train <- data_train[colname_G]
H_train <- data_train[colname_H]


# Assemble different regressor matrices
#X_train <- data.matrix(cbind(D_train, F_train, G_train, H_train))
X_train <- data.matrix(cbind(D_train, G_train, H_train))
Z_train <- data.matrix(cbind(G_train, H_train))
#Z_train <- data.matrix(cbind(F_train, G_train, H_train))



# Grid of lasso penalties
lambda_grid <- seq(0, 2, by=.025)

# 1. Use Lasso to shrink/select covariates given their association with outcome Y
lasso_one <- glmnet(X_train, y_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
# 2. Use Lasso to shrink/select covariates given their association with tratment D
#lasso_two <- glmnet(GH_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
lasso_two <- glmnet(Z_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

# Extract estimates from both models
beta_hats_lasso_one <- t(matrix(coef(lasso_one, s=lambda_grid)[-1,], nrow=dim(X_train)[2], ncol=length(lambda_grid)))
beta_hats_lasso_two <- t(matrix(coef(lasso_two, s=lambda_grid)[-1,], nrow=dim(Z_train)[2], ncol=length(lambda_grid)))



# 1.
df <- data.frame(cbind(beta_hats_lasso_one, lambda_grid))
names(df) <- colnames_DGH
#names(df) <- colnames_DFGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

lasso_coef_shrink_plot(df=df)


# 2.
df <- data.frame(cbind(beta_hats_lasso_two, lambda_grid))
names(df) <- colnames_GH
#names(df) <- colnames_FGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

lasso_coef_shrink_plot(df=df)










#TODO Quick inspection
#plot(lasso_one_cv)
#plot(lasso_one, label=TRUE)
#print(lasso_one)
#coef(lasso_one, s=.1)

# >>>
# Quick look

data_set <- generate_data(n=n,
                          n_F_attr=n_F_attr,
                          n_G_attr=n_G_attr,
                          n_H_attr=n_H_attr,
                          treatment_effect=treatment_effect,
                          unconfoundedness_rate=unconfoundedness_rate)

data_train <- data_set$data
y_train <- data_train$y
D_train <- data_train$D
F_train <- data.matrix(data_train[colname_F])
G_train <- data.matrix(data_train[colname_G])
H_train <- data.matrix(data_train[colname_H])


# Assemble different regressor matrices
#X_train <- data.matrix(cbind(D_train, F_train, G_train, H_train))
X_train <- data.matrix(cbind(D_train, G_train, H_train))
Z_train <- data.matrix(cbind(G_train, H_train))
#Z_train <- data.matrix(cbind(F_train, G_train, H_train))

# Regression results
#lm(y_train ~ D_train + F_train + G_train + H_train)
lm(y_train ~ D_train + G_train + H_train)
lm(y_train ~ D_train + H_train)
lm(y_train ~ D_train)

#lm(D_train ~ F_train + G_train + H_train)
#lm(D_train ~ G_train + H_train)
#lm(D_train ~ G_train)

#beta_hat <- least_squares_estimator(y=y_train, X=X_train)
# <<<




















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
