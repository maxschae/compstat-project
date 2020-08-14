### CompStat project --- Max Schäfer

rm(list=ls())


#install.packages("glmnet")
#install.packages("randomcoloR")

library(glmnet)
library(mvtnorm)

library(stringr)

library(ggplot2)
library(gridExtra)
library(randomcoloR)

# Own modules
#source("data_generating_process.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/data_generating_process.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/plot_functions.R")


# Plotting section ---


### Main ---
set.seed(1234)

n <- 200

n_F_attr <- 30
n_G_attr <- 30
n_H_attr <- 30

beta_F_size <- 1
beta_H_size <- 1

beta_GD_size <- .5
beta_GY_size <- .25
treatment_effect <- .25

unconfoundedness_rate <- .5


# Belloni et al. (2014) setup
n <- 400

n_G_attr <- 20
beta_DY <- .5




# Housekeeping of feature names and identifiers for identification.
# Create vector with covariate names.
colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))

# Required for first selection step
colnames_DG <- c("D", colname_G)
# Required for second selection step
colnames_GH <- c(colname_G, colname_H)
colnumbers_GH <- seq(from=1, to=(n_G_attr+n_H_attr), by=1)


# Required for first selection step
colnames_DG <- c("D", colname_G)
colnames_DGH <- c("D", colname_G, colname_H)

# Required for second selection step
colnames_GH <- c(colname_G, colname_H)


#TODO Do this for multiple runs... quite tricky since its 3-dimensions
# i.e. lambdas, coefs and simulation runs span three dimensions


data_set <- generate_data(n=n,
                          n_F_attr=n_F_attr,
                          n_G_attr=n_G_attr,
                          n_H_attr=n_H_attr,
                          treatment_effect=treatment_effect,
                          beta_GD_size=beta_GD_size,
                          beta_GY_size=beta_GY_size,
                          beta_H_size=beta_H_size,
                          beta_F_size=beta_F_size,
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
#X_train <- data.matrix(cbind(D_train, G_train, H_train))
X_train <- data.matrix(cbind(D_train, G_train))
#Z_train <- data.matrix(cbind(G_train, H_train))
Z_train <- data.matrix(cbind(G_train))
#Z_train <- data.matrix(cbind(F_train, G_train, H_train))


# Get true covariate identifier
true_covariate_identifier <- data_set$true_covariate_identifier
names(true_covariate_identifier) <- colnames_GH
covariate_identifier_G <- true_covariate_identifier[colname_G]
covariate_identifier_G[covariate_identifier_G > 0] <- "n"
covariate_identifier_G[covariate_identifier_G != "n"] <- "p"
# Change column name to identify sparse covariates
colname_G_effect <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1), covariate_identifier_G)
colnames_DG_effect <- c("D", colname_G_effect)
colnames_DGH_effect <- c("D", colname_G_effect, colname_H)
colnames_GH_effect <- c(colname_G_effect, colname_H)



# Grid of lasso penalties
lambda_grid <- seq(0, 2, by=.025)
#lambda_grid <- seq(0, .15, by=.001)

# 1. Use Lasso to shrink/select covariates given their association with outcome Y
lasso_one <- glmnet(X_train, y_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

#TODO
# 1. Use Adaptive Lasso to shrink/select covariates given their association with outcome Y
#beta_hat_ols <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
# Without intercept
#beta_hat_ols <- beta_hat_ols[2:length(beta_hat_ols)]
#beta_hat_ols[1] <- 10
#ad_lasso_one <- glmnet(X_train, y_train, alpha=1, penalty.factor=1/abs(beta_hat_ols),
#                       lambda=lambda_grid, intercept=FALSE)

# 2. Use Lasso to shrink/select covariates given their association with tratment D
#lasso_two <- glmnet(GH_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
lasso_two <- glmnet(Z_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

# Extract estimates from both models
beta_hats_lasso_one <- t(matrix(coef(lasso_one, s=lambda_grid)[-1,], nrow=dim(X_train)[2], ncol=length(lambda_grid)))
#beta_hats_ad_lasso_one <- t(matrix(coef(ad_lasso_one, s=lambda_grid)[-1,], nrow=dim(X_train)[2], ncol=length(lambda_grid)))
beta_hats_lasso_two <- t(matrix(coef(lasso_two, s=lambda_grid)[-1,], nrow=dim(Z_train)[2], ncol=length(lambda_grid)))


#
lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE,
                          type.measure="mse", nfolds=10)
# MSE-optimal penalty term for lasso
lambda_min_one <- lasso_one_cv$lambda.min
lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.min")[-1])

#ad_lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, penalty.factor=1/abs(beta_hat_ols), intercept=FALSE,
#                          type.measure="mse", nfolds=10)
# MSE-optimal penalty term for lasso
#lambda_min_one_ad <- ad_lasso_one_cv$lambda.min

lasso_two_cv <- cv.glmnet(Z_train, D_train, alpha=1, intercept=FALSE,
                          type.measure="mse", nfolds=10)
# MSE-optimal penalty term for lasso
lambda_min_two <- lasso_two_cv$lambda.min
lasso_coefs_two <- cbind(coef(lasso_two_cv, s="lambda.min")[-1])




# 1.
df <- data.frame(cbind(beta_hats_lasso_one, lambda_grid))
#df <- data.frame(cbind(beta_hats_ad_lasso_one, lambda_grid))
names(df) <- colnames_DG_effect
#names(df) <- colnames_DGH_effect
#names(df) <- colnames_DFGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

p1 <- lasso_coef_shrink_plot(df=df, lambda_min=lambda_min_one)
p1



# 2.
df <- data.frame(cbind(beta_hats_lasso_two, lambda_grid))
#names(df) <- colname_G
names(df) <- colname_G_effect
#names(df) <- colnames_GH_effect
#names(df) <- colnames_FGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

p2 <- lasso_coef_shrink_plot(df=df, lambda_min=lambda_min_two)
p2


grid.arrange(p1, p2, nrow=1)



length(lasso_coefs_one[lasso_coefs_one==0])
length(lasso_coefs_two[lasso_coefs_two==0])
unconfoundedness_rate*n_G_attr


#
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/data_generating_process.R")

data_set <- generate_data(n=n,
                          n_F_attr=n_F_attr,
                          n_G_attr=n_G_attr,
                          n_H_attr=n_H_attr,
                          treatment_effect=treatment_effect,
                          beta_GD_size=beta_GD_size,
                          beta_GY_size=beta_GY_size,
                          beta_H_size=beta_H_size,
                          beta_F_size=beta_F_size,
                          unconfoundedness_rate=unconfoundedness_rate)

data <- data_set$data
data_train <- data[1:dim(data)[1]/2, ]

# Draw training data
y_train <- data_train$y
D_train <- data_train$D
F_train <- data_train[colname_F]
G_train <- data_train[colname_G]
H_train <- data_train[colname_H]

G_train <- data.matrix(G_train)
H_train <- data.matrix(H_train)
summary(lm(D_train ~ G_train))$r.squared
summary(lm(y_train ~ D_train + G_train))$r.squared
summary(lm(y_train ~ D_train + G_train + H_train))$r.squared
