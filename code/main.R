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


n <- 200
n_F_attr <- 10
n_G_attr <- 20
n_H_attr <- 30


# Housekeeping of feature names and identifiers for identification.
# Create vector with covariate names.
colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))

# Required for first selection step
colnames_DGH <- c("D", colname_G, colname_H)
colnames_DFGH <- c("D", colname_F, colname_G, colname_H)

# Required for second selection step
colnames_GH <- c(colname_G, colname_H)
colnames_FGH <- c(colname_F, colname_G, colname_H)

covariate_identifier <- seq(from=1, to=(n_F_attr+n_G_attr+n_H_attr), by=1)




### Double-selection method

###TODO write into loop ; or use 'replicate' to make vectorized operation
#, the latter seems not too obvious ...

double_select_covariates <- function(data_train) {

  y_train <- data_train$y
  D_train <- data_train$D

  # Assemble different regressor matrices
  X_train <- cbind(data_train$D, data_train$F, data_train$G, data_train$H)
  #X_train <- cbind(data_train$D, data_train$G, data_train$H)
  GH_train <- cbind(data_train$G, data_train$H)
  FGH_train <- cbind(data_train$F, data_train$G, data_train$H)


  # Double-selection method

  # 1. Use Lasso to select covariates given their association with outcome Y
  lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.min")[-1])
  # Exclude coef. of treatment since we are interested in selecting covariates
  lasso_coefs_one <- lasso_coefs_one[2:dim(X_train)[2], ]

  # Create vector that indicates selected covariates by lasso_one_cv
  selected_covars_one <- seq(from=1, to=length(lasso_coefs_one), by=1)
  selected_covars_one <- t(cbind(selected_covars_one))
  # Name covariates.
  #names(selected_covars_one) <- colnames_GH
  names(selected_covars_one) <- colnames_FGH
  # Keep indicators of selected covariates.
  selected_covars_one[lasso_coefs_one==0] <- NaN
  #selected_covars_one <- selected_covars_one[!is.na(selected_covars_one)]



  # 2. Use Lasso to select covariates given their association with tratment D
  lasso_two_cv <- cv.glmnet(FGH_train, D_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lasso_coefs_two <- cbind(coef(lasso_two_cv, s="lambda.min")[-1])

  # Create vector that indicates selected covariates by lasso_two_cv
  selected_covars_two <- seq(from=1, to=length(lasso_coefs_two), by=1)
  selected_covars_two <- t(cbind(selected_covars_two))
  # Name covariates
  #names(selected_covars_two) <- colnames_GH
  names(selected_covars_two) <- colnames_FGH
  # Keep indicators of selected covariates.
  selected_covars_two[lasso_coefs_two==0] <- NaN
  #selected_covars_two <- selected_covars_two[!is.na(selected_covars_two)]

  #return(list(selection_one=lasso_coefs_one, selection_two=lasso_coefs_two))
  return(list(selection_one=selected_covars_one, selection_two=selected_covars_two))
}



###
R <- 20


for (r in 1:R) {

  # Draw training data
  data_train <- generate_data(n=n, n_F_attr=n_F_attr, n_G_attr=n_G_attr, n_H_attr=n_H_attr)

  # Select covariates according to the double-selection method
  selected_covars <- double_select_covariates(data_train=data_train)
  # Exclude covariates only if both lasso_one and lasso_two suggest it
  double_selected_covars <- rbind(selected_covars$selection_one,
                                  selected_covars$selection_two,
                                  covariate_identifier)

  double_selected_covars <- double_selected_covars[ , colSums(is.na(double_selected_covars)) < 2]
  double_selected_covars <- double_selected_covars[3, ]


  # Draw test data and exclude covariates that were dropped by the double-selection method
  data_test <- generate_data(n=n, n_F_attr=n_F_attr, n_G_attr=n_G_attr, n_H_attr=n_H_attr)

  y_test <- data_test$y
  D_test <- data_test$D

  # Assemble different regressor matrices
  X_test <- cbind(D_test, data_test$F, data_test$G, data_test$H)
  #X_test <- cbind(D_test, G_test, H_test)
  GH_test <- cbind(data_test$G, data_test$H)
  FGH_test <- cbind(data_test$F, data_test$G, data_test$H)



  #TODO evaluate by some metric ... number of falsely dropped covariates


  FGH_test_final <- FGH_test[1:n, double_selected_covars]
  X_test_final <- cbind(D_test, FGH_test_final)



  #TODO evaluate by some metric ... bias of D_coef (causality), MSE (prediction)
  lm_result <- lm(y_test ~ X_test_final)

}


# Collect lasso beta_estimates
#df_selected_covars_one[r, ] <- selected_covars$selection_one
#df_selected_covars_two[r, ] <- selected_covars$selection_two

# OLS with union of X_select_one and X_select_two


#TODO #TODO #TODO #TODO #TODO
# Perhaps, draw dataset and then split into training and test.
# Why? Since otherwise G7 which was dropped based on its association
# in this particular draw, may not be droppable next draw (test_draw)
# or am I missing the whole point about train and test split





# 2. Use Lasso to select covariates given their association with treatment D
#TODO cannot use logistic model then?? or use extension of lasso with logit...
# type.measure class gives misclassification error
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


# Draw training data
y_train <- data_train$y
D_train <- data_train$D

# Assemble different regressor matrices
X_train <- cbind(data_train$D, data_train$F, data_train$G, data_train$H)
#X_train <- cbind(data_train$D, data_train$G, data_train$H)
GH_train <- cbind(data_train$G, data_train$H)
FGH_train <- cbind(data_train$F, data_train$G, data_train$H)



# Grid of lasso penalties
lambda_grid <- seq(0, 2, by=.025)

# 1. Use Lasso to shrink/select covariates given their association with outcome Y
lasso_one <- glmnet(X_train, y_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
# 2. Use Lasso to shrink/select covariates given their association with tratment D
#lasso_two <- glmnet(GH_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
lasso_two <- glmnet(FGH_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

# Extract estimates from both models
beta_hats_lasso_one <- t(matrix(coef(lasso_one, s=lambda_grid)[-1,], nrow=dim(X_train)[2], ncol=length(lambda_grid)))
beta_hats_lasso_two <- t(matrix(coef(lasso_two, s=lambda_grid)[-1,], nrow=dim(FGH_train)[2], ncol=length(lambda_grid)))



# 1.
df <- data.frame(cbind(beta_hats_lasso_one, lambda_grid))
#names(df) <- colnames_DGH
names(df) <- colnames_DFGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

lasso_coef_shrink_plot(df=df)


# 2.
df <- data.frame(cbind(beta_hats_lasso_two, lambda_grid))
#names(df) <- colnames_GH
names(df) <- colnames_FGH
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
data <- generate_data(n=n, n_F_attr=n_F_attr, n_G_attr=n_G_attr, n_H_attr=n_H_attr)
y <- data$y
D <- data$D
F <- data$F
G <- data$G
H <- data$H

# Regression results
lm(y ~ D + F + G + H)
#lm(y ~ D + G + H)
#lm(y ~ D + H)
lm(y ~ D)

lm(D ~ F + G + H)
#lm(D ~ G + H)
lm(D ~ G)

#beta_hat <- least_squares_estimator(y=y, X=cbind(D, F, G, H))
#beta_hat <- least_squares_estimator(y=y, X=cbind(D, G, H))
#beta_hat <- least_squares_estimator(y=y, X=cbind(D, F, H))
#beta_hat <- least_squares_estimator(y=y, X=cbind(D, H))
#beta_hat <- least_squares_estimator(y=y, X=D)
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
