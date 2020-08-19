### CompStat project --- Max Schäfer
# This script contains the variable selection methods used


simple_select_covariates <- function(y_train, X_train) {
  # Simple selection method
  # Off-the-shelve MSE-optimal Lasso; outcome regressed on entire feature matrix

  # Use Lasso to select covariates given their association with outcome Y
  lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)

  lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.min")[-1])
  lasso_coefs_one <- cbind(coef(lasso_one_cv, s="lambda.1se")[-1])

  # Exclude coef. of treatment since we are interested in selecting covariates
  lasso_coefs_one <- lasso_coefs_one[2:dim(X_train)[2], ]

  # Create vector that indicates selected covariates by lasso_one_cv
  selected_covars_one <- rbind(seq(from=1, to=length(lasso_coefs_one), by=1))
  # Keep indicators of selected covariates.
  selected_covars_one[lasso_coefs_one==0] <- NaN

  return(selected_covars_one)
}



double_select_covariates <- function(D_train, Z_train, y_train, X_train) {
  # Double-selection procedure following Belloni et al. (2014)

  # 1. Use Lasso to select covariates given their association with outcome Y
  selected_covars_one <- simple_select_covariates(y_train=y_train, X_train=X_train)

  # 2. Use Lasso to select covariates given their association with tratment D
  lasso_two_cv <- cv.glmnet(Z_train, D_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lasso_coefs_two <- cbind(coef(lasso_two_cv, s="lambda.min")[-1])
  lasso_coefs_two <- cbind(coef(lasso_two_cv, s="lambda.1se")[-1])

  # Create vector that indicates selected covariates by lasso_two_cv
  selected_covars_two <- rbind(seq(from=1, to=length(lasso_coefs_two), by=1))
  # Keep indicators of selected covariates.
  selected_covars_two[lasso_coefs_two==0] <- NaN

  # Exclude covariates only if both lasso_one and lasso_two suggest it
  colnumbers <- seq(from=1, to=dim(Z_train)[2], by=1)
  double_selected_identifier <- rbind(selected_covars_one,
                                      selected_covars_two,
                                      colnumbers)

  double_selected_identifier[3, colSums(is.na(double_selected_identifier)) == 2] <- NaN
  double_selected_identifier <- double_selected_identifier[3, ]

  return(double_selected_identifier)
}
