### CompStat project --- Max Schäfer
# This script contains the variable selection methods used


simple_select_covariates <- function(data_train) {
  # Simple selection method
  # Off-the-shelve MSE-optimal Lasso; outcome regressed on entire feature matrix

  y_train <- data_train$y
  D_train <- data_train$D
  F_train <- data_train[colname_F]
  G_train <- data_train[colname_G]
  H_train <- data_train[colname_H]

  # Assemble regressor matrix
  X_train <- data.matrix(cbind(D_train, G_train, H_train))
  X_train <- data.matrix(cbind(D_train, G_train))

  # Use Lasso to select covariates given their association with outcome Y
  lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)

  lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.min")[-1])
  lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.1se")[-1])

  select_treatment_effect <- 0 #TODO
  if (lasso_coefs_one[1, ] == 0) {
    select_treatment_effect <- 1
  }

  # Exclude coef. of treatment since we are interested in selecting covariates
  lasso_coefs_one <- lasso_coefs_one[2:dim(X_train)[2], ]

  # Create vector that indicates selected covariates by lasso_one_cv
  selected_covars_one <- rbind(seq(from=1, to=length(lasso_coefs_one), by=1))
  # Name covariates for identification later
  #names(selected_covars_one) <- colnames_GH
  names(selected_covars_one) <- colname_G
  # Keep indicators of selected covariates.
  selected_covars_one[lasso_coefs_one==0] <- NaN

  #return(selected_covars_one)
  return(list(selected_covars_one=selected_covars_one, select_treatment_effect=select_treatment_effect))
}



double_select_covariates <- function(data_train) {
  # Double-selection procedure following Belloni et al. (2014)

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

  X_train <- data.matrix(cbind(D_train, G_train))
  Z_train <- data.matrix(G_train)


  # Double-selection method

  # 1. Use Lasso to select covariates given their association with outcome Y
  selected_covars_one <- simple_select_covariates(data_train=data_train)$selected_covars_one

  # 2. Use Lasso to select covariates given their association with tratment D
  lasso_two_cv <- cv.glmnet(Z_train, D_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lasso_coefs_two <- cbind(coef(lasso_two_cv, s="lambda.min")[-1])
  lasso_coefs_two <- cbind(coef(lasso_two_cv, s="lambda.1se")[-1])

  # Create vector that indicates selected covariates by lasso_two_cv
  selected_covars_two <- rbind(seq(from=1, to=length(lasso_coefs_two), by=1))
  # Name covariates
  #names(selected_covars_two) <- colnames_GH
  names(selected_covars_two) <- colname_G
  # Keep indicators of selected covariates.
  selected_covars_two[lasso_coefs_two==0] <- NaN

  # Exclude covariates only if both lasso_one and lasso_two suggest it
  double_selected_identifier <- rbind(selected_covars_one,
                                      selected_covars_two,
                                      colnumbers_G)
                                      #colnumbers_GH

  double_selected_identifier[3, colSums(is.na(double_selected_identifier)) == 2] <- NaN
  double_selected_identifier <- double_selected_identifier[3, ]

  return(double_selected_identifier)
}
