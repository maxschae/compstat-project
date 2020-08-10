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
  select_treatment_effect <- 0 #TODO
  if (lasso_coefs_one[1, ] == 0) {
    select_treatment_effect <- 1
  }
  #TODO

  lasso_coefs_one <- lasso_coefs_one[2:dim(X_train)[2], ]

  # Create vector that indicates selected covariates by lasso_one_cv
  selected_covars_one <- seq(from=1, to=length(lasso_coefs_one), by=1)
  selected_covars_one <- rbind(selected_covars_one)
  # Name covariates.
  names(selected_covars_one) <- colnames_GH
  #names(selected_covars_one) <- colnames_FGH
  # Keep indicators of selected covariates.
  selected_covars_one[lasso_coefs_one==0] <- NaN

  #return(selected_covars_one)
  return(list(selected_covars_one=selected_covars_one, select_treatment_effect=select_treatment_effect))
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
  selected_covars_one <- simple_select_covariates(data_train=data_train)$selected_covars_one

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
                                    beta_GD_size=1, beta_GY_size=1, beta_H_size=1,
                                    unconfoundedness_rate, selection_method) {


  # Draw data once and split into train and test (*)
  data_set <- generate_data(n=n,
                        n_F_attr=n_F_attr,
                        n_G_attr=n_G_attr,
                        n_H_attr=n_H_attr,
                        treatment_effect=treatment_effect,
                        beta_GD_size=beta_GD_size,
                        beta_GY_size=beta_GY_size,
                        beta_H_size=beta_H_size,
                        unconfoundedness_rate=unconfoundedness_rate)

  data <- data_set$data

  # Shuffle data before train-test-split to ensure no unintended structure exists
  # due to the order of observations
  random_rows <- sample(nrow(data))
  data <- data[random_rows, ]
  # Split data into training and test sets
  data_train <- data[1:(dim(data)[1]/2), ]
  data_test <- data[(dim(data)[1]/2+1):dim(data)[1], ]

  true_covariate_identifier <- data_set$true_covariate_identifier

  # *Selection section*
  # Select covariates according to the simple or double-selection method
  if (selection_method=="simple") {
    #selection_identifier <- simple_select_covariates(data_train=data_train)
    out <- simple_select_covariates(data_train=data_train)
    selection_identifier <- out$selected_covars_one
    select_treatment_identifier <- out$select_treatment_effect
  }
  if (selection_method=="double") {
    selection_identifier <- double_select_covariates(data_train=data_train)
  }


  covariate_identifier <- data.frame(rbind(true_covariate_identifier,
                                           selection_identifier))

  #names(covariate_identifier) <- colnames_FGH
  names(covariate_identifier) <- colnames_GH

  #covariate_identifier_F <- covariate_identifier[colname_F] #TODO F_
  covariate_identifier_G <- covariate_identifier[colname_G]
  covariate_identifier_H <- covariate_identifier[colname_H]


  #TODO check whether treatment was also dropped?

  # The setup here is labeled as follows:
  # . 'positive' refers to state where covariate has a zero effect.
  # . 'negative' refers to state where covariate has a non-zero effect.
  # Thus, true-positive refers to state where covariate has a zero effect
  # and the selection methods excludes covariate.
  # True-negative refers to state where covariate has a non-zero effect
  # and the selection method does not exclude covariate.
  # False-positive : covariate has a zero effect and is not excluded.
  # False-negative : covariate has non-zero effect and is not excluded.

  tp_selection_count_G <- 0
  tn_selection_count_G <- 0
  fp_selection_count_G <- 0
  fn_selection_count_G <- 0

  covariate_identifier_G[is.na(covariate_identifier_G)] <- 0

  for (j in 1:dim(covariate_identifier_G)[2]) {
    # True : the correct choice was made
    if (covariate_identifier_G[1, j] == covariate_identifier_G[2, j]) {
      # True positive : an un-associated variable was excluded
      if (covariate_identifier_G[1, j] == 0) {
        tp_selection_count_G <- tp_selection_count_G + 1
      }
      # True negative : an associated variable was not excluded
      if (covariate_identifier_G[1, j] != 0) {
        tn_selection_count_G <- tn_selection_count_G + 1
      }
    }
    # False : the false choice was made
    if (covariate_identifier_G[1, j] != covariate_identifier_G[2, j]) {
      # False positive : an associated variable was excluded
      if (covariate_identifier_G[1, j] != 0) {
        fp_selection_count_G <- fp_selection_count_G + 1
      }
      # False negative : an un-associated variable was not excluded
      if (covariate_identifier_G[1, j] == 0) {
        fn_selection_count_G <- fn_selection_count_G + 1
      }
    }
  }


  #TODO
  tp_selection_rate_G <- tp_selection_count_G #/ n_G_attr #/(unconfoundedness_rate+1e-4)
  tn_selection_rate_G <- tn_selection_count_G #/ n_G_attr #/(unconfoundedness_rate+1e-4)
  fp_selection_rate_G <- fp_selection_count_G #/ n_G_attr #/(unconfoundedness_rate+1e-4)
  fn_selection_rate_G <- fn_selection_count_G #/ n_G_attr #/(unconfoundedness_rate+1e-4)




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

  Z_test_selected <- Z_test[, selection_identifier[!is.na(selection_identifier)]]
  X_test_selected <- data.matrix(cbind(D_test, Z_test_selected))

  # Estimate treatment effect by standard OLS provided with double-selected covariates
  lm_result <- lm(y_test ~ X_test_selected)
  # Collect metrics of interest
  confounder_bias <- lm_result$coefficients[2] - treatment_effect
  test_MSE <- (1/length(lm_result$residuals)) * sum((lm_result$residuals)**2)

  return(list(confounder_bias=confounder_bias,
              test_MSE=test_MSE,
              tp_selection_rate_G=tp_selection_rate_G,
              tn_selection_rate_G=tn_selection_rate_G,
              fp_selection_rate_G=fp_selection_rate_G,
              fn_selection_rate_G=fn_selection_rate_G,
              select_treatment_identifier=select_treatment_identifier))
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




### Main ---
set.seed(1234)


n <- 400
n_F_attr <- 70
n_G_attr <- 90
n_H_attr <- 21
treatment_effect <- 1
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



# Main : run simulation
R <- 20

sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                       n_F_attr=n_F_attr,
                                       n_G_attr=n_G_attr,
                                       n_H_attr=n_H_attr,
                                       treatment_effect=treatment_effect,
                                       unconfoundedness_rate=unconfoundedness_rate,
                                       selection_method="simple"))

confounder_bias_vec <- unlist(sim_results_vec[1, 1:R])
test_MSE_vec <- unlist(sim_results_vec[2, 1:R])

tp_selection_rate_G <- unlist(sim_results_vec[3, 1:R])
tn_selection_rate_G <- unlist(sim_results_vec[4, 1:R])
fp_selection_rate_G <- unlist(sim_results_vec[5, 1:R])
fn_selection_rate_G <- unlist(sim_results_vec[6, 1:R])

select_treatment_identifier <- unlist(sim_results_vec[7, 1:R]) #TODO

# Evaluate
cat("Mean test MSE is", mean(test_MSE_vec))
cat("Median confounder bias is", median(confounder_bias_vec))

hist(confounder_bias_vec)




# Simulations
#TODO iterating over treatment_effect or numbers of attributes or unconfoundedness rate / intensity etc.
# this forms the key simulation study ---
# Setup
selection_method <- "simple" #double
unconfoundedness_rate <- .5
treatment_effect <- 1

# Vary over
iter_over <- "n_G_attr"

unconfoundedness_rate_vec <- seq(from=0, to=1, by=.1)
treatment_effect_vec <- seq(from=0, to=-3, by=-.25)
beta_GD_size_vec <- seq(from=0, to=1.5, by=.05)
beta_GY_size_vec <- seq(from=0, to=2, by=.05)
beta_H_size_vec <- seq(from=0, to=4, by=.1)

n_G_attr_vec <- seq(from=10, to=120, by=10)

if (iter_over=="unconfoundedness_rate") {
  sim_parameter_vec <- unconfoundedness_rate_vec
}
if (iter_over=="treatment_effect") {
  sim_parameter_vec <- treatment_effect_vec
}
if (iter_over=="n_G_attr") {
  sim_parameter_vec <- n_G_attr_vec
}

#TODO
#sim_parameter_vec <- beta_GY_size_vec
#sim_parameter_vec <- beta_GD_size_vec
#sim_parameter_vec <- beta_H_size_vec

#n_H_attr <- 30  #TODO
#n_G_attr <- 30  #TODO
#n_F_attr <- 30  #TODO



#print("Simulation parameter varied over is: ")

# Initialize containers
confounder_bias_mat <- matrix(NA, nrow=length(sim_parameter_vec), ncol=R)
names(confounder_bias_mat) <- str_c(rep("run_", R), seq(from=1, to=R, by=1))

mean_squared_bias_vec <- rep(NA, length(sim_parameter_vec))
mean_abs_dev_bias_vec <- rep(NA, length(sim_parameter_vec))
median_bias_vec <- rep(NA, length(sim_parameter_vec))
mean_test_MSE_vec <- rep(NA, length(sim_parameter_vec))
tp_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))
tn_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))
fp_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))
fn_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))

select_treatment_identifier_vec <- rep(NA, length(sim_parameter_vec))



i <- 1

for (SIM_PARAMETER in sim_parameter_vec) {

  if (iter_over=="unconfoundedness_rate") {
    sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                         n_F_attr=n_F_attr,
                                         n_G_attr=n_G_attr,
                                         n_H_attr=n_H_attr,
                                         treatment_effect=treatment_effect,
                                         unconfoundedness_rate=SIM_PARAMETER,
                                         selection_method=selection_method))
  }
  if (iter_over=="treatment_effect") {
    sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                         n_F_attr=n_F_attr,
                                         n_G_attr=n_G_attr,
                                         n_H_attr=n_H_attr,
                                         treatment_effect=SIM_PARAMETER,
                                         unconfoundedness_rate=unconfoundedness_rate,
                                         selection_method=selection_method))
  }
  if (iter_over=="beta_GY_size") {
    sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                         n_F_attr=n_F_attr,
                                         n_G_attr=n_G_attr,
                                         n_H_attr=n_H_attr,
                                         treatment_effect=treatment_effect,
                                         beta_GY_size=SIM_PARAMETER,
                                         unconfoundedness_rate=unconfoundedness_rate,
                                         selection_method=selection_method))
  }
  if (iter_over=="beta_GD_size") {
    sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                         n_F_attr=n_F_attr,
                                         n_G_attr=n_G_attr,
                                         n_H_attr=n_H_attr,
                                         treatment_effect=treatment_effect,
                                         beta_GD_size=SIM_PARAMETER,
                                         unconfoundedness_rate=unconfoundedness_rate,
                                         selection_method=selection_method))
  }
  if (iter_over=="beta_H_size") {
    sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                         n_F_attr=n_F_attr,
                                         n_G_attr=n_G_attr,
                                         n_H_attr=n_H_attr,
                                         treatment_effect=treatment_effect,
                                         beta_H_size=SIM_PARAMETER,
                                         unconfoundedness_rate=unconfoundedness_rate,
                                         selection_method=selection_method))
  }
  if (iter_over=="n_G_attr") {
    colname_G <- str_c(rep("G_", SIM_PARAMETER), seq(from=1, to=SIM_PARAMETER, by=1))
    colnames_GH <- c(colname_G, colname_H)
    sim_results_vec <- replicate(R, compute_confounder_bias(n=n,
                                         n_F_attr=n_F_attr,
                                         n_G_attr=SIM_PARAMETER,
                                         n_H_attr=n_H_attr,
                                         treatment_effect=treatment_effect,
                                         unconfoundedness_rate=unconfoundedness_rate,
                                         selection_method=selection_method))
  }

  # Metrics --- confounder bias
  confounder_bias_mat[i, ] <- unlist(sim_results_vec[1, 1:R])
  mean_squared_bias_vec[i] <- mean((unlist(sim_results_vec[1, 1:R]))**2)
  mean_abs_dev_bias_vec[i] <- mean(abs(unlist(sim_results_vec[1, 1:R])))
  median_bias_vec[i] <- median(unlist(sim_results_vec[1, 1:R]))
  # Metrics --- prediction
  mean_test_MSE_vec[i] <- mean(unlist(sim_results_vec[2, 1:R]))
  # Metrics --- variable selection
  tp_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[3, 1:R]))
  tn_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[4, 1:R]))
  fp_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[5, 1:R]))
  fn_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[6, 1:R]))

  select_treatment_identifier_vec[i] <- mean(unlist(sim_results_vec[7, 1:R]))

  i <- i+1
}

# If zero, D was *not* excluded
select_treatment_identifier_vec

# Evaluation
median_bias_vec <- apply(confounder_bias_mat, MARGIN=1, FUN=median)




### Visualization
# Metrics --- confounder bias
x_axis <- sim_parameter_vec
df <- data.frame(cbind(x_axis, mean_squared_bias_vec, mean_abs_dev_bias_vec, median_bias_vec))

# Plot
plot <- ggplot(data=df, aes(x=x_axis))
plot <- plot +
        geom_line(aes_string(y="mean_squared_bias_vec"), size=1, col="blue", alpha=1) +
        geom_line(aes_string(y="mean_abs_dev_bias_vec"), size=1, col="red", alpha=1) +
        geom_line(aes_string(y="median_bias_vec"), size=1, col="grey", alpha=1) +
        geom_line(aes(y=0), lty="dashed") +
        ylab("Bias of Treatment Coefficient") +
        xlab("n_G_attr")

plot


# Metrics --- prediction
x_axis <- sim_parameter_vec
df <- data.frame(cbind(x_axis, mean_test_MSE_vec))

# Plot
plot <- ggplot(data=df, aes(x=x_axis))
plot <- plot +
        geom_line(aes_string(y="mean_test_MSE_vec"), size=1, col="blue", alpha=1) +
        geom_line(aes(y=0), lty="dashed") +
        ylab("Mean test MSE (prediction)") +
        xlab("n_G_attr")

plot


# Metrics --- variable selection
x_axis <- sim_parameter_vec
df <- data.frame(cbind(x_axis,
                       tp_selection_rate_G_vec, tn_selection_rate_G_vec,
                       fp_selection_rate_G_vec, fn_selection_rate_G_vec))

# Plot
plot <- ggplot(data=df, aes(x=x_axis))
plot <- plot +
        geom_line(aes_string(y="tp_selection_rate_G_vec"), size=1, col="blue", alpha=.9) +
        geom_line(aes_string(y="tn_selection_rate_G_vec"), size=1, col="blue", alpha=.5) +
        geom_line(aes_string(y="fp_selection_rate_G_vec"), size=1, col="red", alpha=.9) +
        geom_line(aes_string(y="fn_selection_rate_G_vec"), size=1, col="red", alpha=.5) +
        geom_line(aes(y=0), lty="dashed") +
        ylab("TP, TN, FP, FN rates for selecting confounder") +
        xlab("n_G_attr")

plot




library(tidyr)

df_long <- tidyr::gather(df[ , 2:dim(df)[2]], key="metric")
df_long["x_axis"] <- df[,1]
plot <- ggplot(data=df_long, aes(x=x_axis, y=value, colour=metric)) +
        geom_line() +
        geom_line(aes(y=0), lty="dashed") +
        ylab("TP, TN, FP, FN rates for selecting confounder") +
        xlab("Unconfoundedness Rate")

plot



# 2. Use Lasso to select covariates given their association with treatment D
#TODO cannot use logistic model then?? or use extension of lasso with logit...
# type.measure class gives misclassification error for cv
#lasso_two_cv <- cv.glmnet(G_train, D_train, alpha=1, family="binomial", type.measure="class", intercept=FALSE)














#TODO Quick inspection
#plot(lasso_one_cv)
#plot(lasso_one, label=TRUE)
#print(lasso_one)
#coef(lasso_one, s=.1)

# >>>
# Quick look
n_F_attr <- 12
n_G_attr <- 12
n_H_attr <- 12
treatment_effect <- 1
unconfoundedness_rate <- 0


colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))


data_set <- generate_data(n=n,
                          n_F_attr=n_F_attr,
                          n_G_attr=n_G_attr,
                          n_H_attr=n_H_attr,
                          treatment_effect=treatment_effect,
                          beta_GY_size = 2,
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
lm(y_train ~ G_train + H_train)
lm(y_train ~ D_train + H_train)

lm(y_train ~ D_train + G_train[, 1:12] + H_train)
lm(y_train ~ D_train)

#lm(D_train ~ F_train + G_train + H_train)
lm(D_train ~ G_train + H_train)
lm(D_train ~ G_train)

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
