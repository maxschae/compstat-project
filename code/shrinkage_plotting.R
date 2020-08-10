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
#source("data_generating_process.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/data_generating_process.R")



# Plot all lasso coefficients with increasing lambda
lasso_coef_shrink_plot <- function(df, lambda_min=NaN) {

  cols <- names(df)[1:(length(names(df))-1)]

  # Initialize plot.
  plot <- ggplot(data=df, aes(x=lambda_grid))
  # Loop through each beta and plot its value with increasing penalty.
  for (name in cols) {
    random_transparency <- runif(n=1, min=.5, max=1) #

    if (startsWith(name, "F")) {
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
    else if (startsWith(name, "D")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#43f7d0")
    }
  }
  plot <- plot +
          geom_line(aes(y=0), lty="dashed") +
          ylab("Lasso coefficients") +
          xlab("Lambda")

  if (is.na(lambda_min) == FALSE) {
    plot <- plot + geom_vline(xintercept=lambda_min, linetype="dashed",
                              color="black", alpha=.8, size=1)
  }
  plot
}


# Plotting section ---


### Main ---
set.seed(1234)

n <- 400
n_F_attr <- 70
n_G_attr <- 60
n_H_attr <- 60
treatment_effect <- -1
unconfoundedness_rate <- .5



#n_F_attr <- 12
n_G_attr <- 120
n_H_attr <- 12
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

#TODO Do this for multiple runs... quite tricky since its 3-dimensions
# i.e. lambdas, coefs and simulation runs span three dimensions


data_set <- generate_data(n=n,
                          n_F_attr=n_F_attr,
                          n_G_attr=n_G_attr,
                          n_H_attr=n_H_attr,
                          treatment_effect=treatment_effect,
                          beta_GD_size=1,
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
lambda_grid <- seq(0, .15, by=.001)

# 1. Use Lasso to shrink/select covariates given their association with outcome Y
lasso_one <- glmnet(X_train, y_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
# 2. Use Lasso to shrink/select covariates given their association with tratment D
#lasso_two <- glmnet(GH_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)
lasso_two <- glmnet(Z_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

# Extract estimates from both models
beta_hats_lasso_one <- t(matrix(coef(lasso_one, s=lambda_grid)[-1,], nrow=dim(X_train)[2], ncol=length(lambda_grid)))
beta_hats_lasso_two <- t(matrix(coef(lasso_two, s=lambda_grid)[-1,], nrow=dim(Z_train)[2], ncol=length(lambda_grid)))


#
lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE,
                          type.measure="mse", nfolds=10)
# MSE-optimal penalty term for lasso
lambda_min_one <- lasso_one_cv$lambda.min
lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.min")[-1])



# 1.
df <- data.frame(cbind(beta_hats_lasso_one, lambda_grid))
names(df) <- colnames_DGH
#names(df) <- colnames_DFGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

lasso_coef_shrink_plot(df=df, lambda_min=lambda_min_one)



# 2.
df <- data.frame(cbind(beta_hats_lasso_two, lambda_grid))
names(df) <- colnames_GH
#names(df) <- colnames_FGH
#TODO Drop treatment for visualization purposes
#df <- subset(df, select=-c(D))

lasso_coef_shrink_plot(df=df)
