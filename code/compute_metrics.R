### CompStat project --- Max Schäfer
# This script contains the function computing all metrics used
# to evaluate the variable selection methods



compute_metrics <- function(data, true_covariate_identifier, selection_method,
                            colnames_covariates, colnames_confounders,
                            treatment_effect) {



  # Shuffle data before train-test-split to ensure no unintended structure exists
  # due to the order of observations
  random_rows <- sample(nrow(data))
  data <- data[random_rows, ]

  # Split data into training and test sets
  data_train <- data[1:(dim(data)[1]/2), ]
  data_test <- data[(dim(data)[1]/2+1):dim(data)[1], ]

  y_train <- data_train$y
  D_train <- data_train$D

  # Assemble regressor matrices for first and second stage
  Z_train <- data.matrix(data_train[colnames_covariates])
  X_train <- data.matrix(cbind(D_train, Z_train))

  # *Selection section*
  # Select covariates according to the SIMPLE selection method
  selection_method <- "simple"
  selection_identifier_simple <- simple_select_covariates(y_train=y_train, X_train=X_train)

  if (length(true_covariate_identifier) != length(selection_identifier_simple)) {
    print("Warning: check how true_covariate_identifier vector is initialized in data_generating_process.R")
  }
  # Collect information which potential confounders were selected
  covariate_identifier <- data.frame(rbind(true_covariate_identifier,
                                           selection_identifier_simple))

  names(covariate_identifier) <- colnames_covariates
  # Select confounders
  covariate_identifier_G <- covariate_identifier[colnames_confounders]
  covariate_identifier_G[is.na(covariate_identifier_G)] <- 0  # work with zero

  # The setup here is labeled as follows:
  # . 'positive' refers to state where covariate has a zero effect.
  # . 'negative' refers to state where covariate has a non-zero effect.
  # Thus, true-positive refers to state where covariate has a zero effect
  # and the selection methods excludes covariate.
  # True-negative refers to state where covariate has a non-zero effect
  # and the selection method does not exclude covariate.
  # False-positive : covariate has a zero effect and is not excluded.
  # False-negative : covariate has non-zero effect and is not excluded.

  # True positive : an un-associated variable was excluded
  tp_selection_count_G <- sum((covariate_identifier_G[1,] == covariate_identifier_G[2,]) & (covariate_identifier_G[1,] == 0))
  # True negative : an associated variable was not excluded
  tn_selection_count_G <- sum((covariate_identifier_G[1,] == covariate_identifier_G[2,]) & (covariate_identifier_G[1,] != 0))
  # False positive : an associated variable was excluded
  fp_selection_count_G <- sum((covariate_identifier_G[1,] != covariate_identifier_G[2,]) & (covariate_identifier_G[1,] != 0))
  # False negative : an un-associated variable was not excluded
  fn_selection_count_G <- sum((covariate_identifier_G[1,] != covariate_identifier_G[2,]) & (covariate_identifier_G[1,] == 0))


  #TODO
  tp_selection_rate_G_simple <- tp_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  tn_selection_rate_G_simple <- tn_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  fp_selection_rate_G_simple <- fp_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  fn_selection_rate_G_simple <- fn_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)


  # Check whether individual confounder was excluded from model
  # For 'houseprices', check first confounder 'purchasing_power'
  #TODO
  selection_confounder_identifier <- selection_identifier_simple
  selection_confounder_identifier[is.na(selection_confounder_identifier)] <- 0
  selection_confounder_identifier_simple <- selection_confounder_identifier[1]


  # Select covariates according to the DOUBLE selection method
  selection_method <- "double"
  selection_identifier_double <- double_select_covariates(D_train=D_train,
                                                          Z_train=Z_train,
                                                          y_train=y_train,
                                                          X_train=X_train)

  if (length(true_covariate_identifier) != length(selection_identifier_double)) {
    print("Warning: check how true_covariate_identifier vector is initialized in data_generating_process.R")
  }
  # Collect information which potential confounders were selected
  covariate_identifier <- data.frame(rbind(true_covariate_identifier,
                                           selection_identifier_double))

  names(covariate_identifier) <- colnames_covariates
  # Select confounders
  covariate_identifier_G <- covariate_identifier[colnames_confounders]
  covariate_identifier_G[is.na(covariate_identifier_G)] <- 0  # work with zero

  tp_selection_count_G <- sum((covariate_identifier_G[1,] == covariate_identifier_G[2,]) & (covariate_identifier_G[1,] == 0))
  tn_selection_count_G <- sum((covariate_identifier_G[1,] == covariate_identifier_G[2,]) & (covariate_identifier_G[1,] != 0))
  fp_selection_count_G <- sum((covariate_identifier_G[1,] != covariate_identifier_G[2,]) & (covariate_identifier_G[1,] != 0))
  fn_selection_count_G <- sum((covariate_identifier_G[1,] != covariate_identifier_G[2,]) & (covariate_identifier_G[1,] == 0))
  #TODO
  tp_selection_rate_G_double <- tp_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  tn_selection_rate_G_double <- tn_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  fp_selection_rate_G_double <- fp_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  fn_selection_rate_G_double <- fn_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)

  # Check whether individual confounder was excluded from model
  # For 'houseprices', check first confounder 'purchasing_power'
  #TODO
  selection_confounder_identifier <- selection_identifier_double
  selection_confounder_identifier[is.na(selection_confounder_identifier)] <- 0
  selection_confounder_identifier_double <- selection_confounder_identifier[1]



  # Compute metrics of interest concerned with bias with test data.
  y_test <- data_test$y
  D_test <- data_test$D
  # Assemble different regressor matrices
  Z_test <- data.matrix(data_test[colnames_covariates])


  # --- Compute metrics for SIMPLE selection
  Z_test_selected_simple <- Z_test[, selection_identifier_simple[!is.na(selection_identifier_simple)]]
  X_test_selected_simple <- data.matrix(cbind(D_test, Z_test_selected_simple))
  # Estimate treatment effect by OLS provided with selected covariates
  lm_result_simple <- lm(y_test ~ X_test_selected_simple)
  # Collect metrics of interest
  confounder_bias_simple <- lm_result_simple$coefficients[2] - treatment_effect
  test_MSE_simple <- (1/length(lm_result_simple$residuals)) * sum((lm_result_simple$residuals)**2)

  # --- Compute metrics for DOUBLE selection
  Z_test_selected_double <- Z_test[, selection_identifier_double[!is.na(selection_identifier_double)]]
  X_test_selected_double <- data.matrix(cbind(D_test, Z_test_selected_double))
  # Estimate treatment effect by OLS provided with selected covariates
  lm_result_double <- lm(y_test ~ X_test_selected_double)
  # Collect metrics of interest
  confounder_bias_double <- lm_result_double$coefficients[2] - treatment_effect
  test_MSE_double <- (1/length(lm_result_double$residuals)) * sum((lm_result_double$residuals)**2)



  return(list(confounder_bias_simple=confounder_bias_simple,
              test_MSE_simple=test_MSE_simple,
              tp_selection_rate_G_simple=tp_selection_rate_G_simple,
              tn_selection_rate_G_simple=tn_selection_rate_G_simple,
              fp_selection_rate_G_simple=fp_selection_rate_G_simple,
              fn_selection_rate_G_simple=fn_selection_rate_G_simple,
              selection_confounder_identifier_simple=selection_confounder_identifier_simple,
              confounder_bias_double=confounder_bias_double,
              test_MSE_double=test_MSE_double,
              tp_selection_rate_G_double=tp_selection_rate_G_double,
              tn_selection_rate_G_double=tn_selection_rate_G_double,
              fp_selection_rate_G_double=fp_selection_rate_G_double,
              fn_selection_rate_G_double=fn_selection_rate_G_double,
              selection_confounder_identifier_double=selection_confounder_identifier_double))
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



compute_metrics_old <- function(dgp, n=200, n_F_attr=30, n_G_attr=30, n_H_attr=30,
                            corr_G=0, treatment_effect=.5,
                            beta_GD_size=.5, beta_GY_size=.25, beta_H_size=.5,
                            beta_F_size=.5, colnames_covariates, colnames_confounders,
                            nonzero_controls=5, selection_method="simple") {


  # Draw data once and split into train and test (*)
  if (dgp == "A") {
    data_set <- generate_data_A(n=n,
                                n_G_attr=n_G_attr,
                                corr_G=corr_G,
                                treatment_effect=treatment_effect,
                                beta_GD_size=beta_GD_size,
                                beta_GY_size=beta_GY_size,
                                nonzero_controls=nonzero_controls)
  }
  if (dgp == "B") {
    data_set <- generate_data_B(n=n,
                                n_G_attr=n_G_attr,
                                n_F_attr=n_F_attr,
                                n_H_attr=n_H_attr,
                                corr_G=corr_G,
                                treatment_effect=treatment_effect,
                                beta_GD_size=beta_GD_size,
                                beta_GY_size=beta_GY_size,
                                beta_H_size=beta_H_size,
                                beta_F_size=beta_F_size,
                                nonzero_controls=nonzero_controls)
  }
  if (dgp == "houseprices") {
    data_set <- generate_data_houseprices_dgp0(n=n, n_G_attr=n_G_attr,
                                               corr_G=corr_G,
                                               beta_GY_inflator=beta_GY_size)
  }

  data <- data_set$data
  # Shuffle data before train-test-split to ensure no unintended structure exists
  # due to the order of observations
  random_rows <- sample(nrow(data))
  data <- data[random_rows, ]

  # Split data into training and test sets
  data_train <- data[1:(dim(data)[1]/2), ]
  data_test <- data[(dim(data)[1]/2+1):dim(data)[1], ]

  y_train <- data_train$y
  D_train <- data_train$D

  # Assemble regressor matrices for first and second stage
  Z_train <- data.matrix(data_train[colnames_covariates])
  X_train <- data.matrix(cbind(D_train, Z_train))

  # *Selection section*
  # Select covariates according to the simple or double-selection method
  if (selection_method=="simple") {
    selection_identifier <- simple_select_covariates(y_train=y_train, X_train=X_train)
  }
  if (selection_method=="double") {
    selection_identifier <- double_select_covariates(D_train=D_train,
                                                     Z_train=Z_train,
                                                     y_train=y_train,
                                                     X_train=X_train)
  }

  # Collect vector identifying non-zero controls (i.e. confounders)
  true_covariate_identifier <- data_set$true_covariate_identifier
  if (length(true_covariate_identifier) != length(selection_identifier)) {
    print("Warning: check how true_covariate_identifier vector is initialized in data_generating_process.R")
  }
  # Collect information which potential confounders were selected
  covariate_identifier <- data.frame(rbind(true_covariate_identifier,
                                           selection_identifier))

  names(covariate_identifier) <- colnames_covariates
  # Select confounders
  covariate_identifier_G <- covariate_identifier[colnames_confounders]
  covariate_identifier_G[is.na(covariate_identifier_G)] <- 0  # work with zero

  # The setup here is labeled as follows:
  # . 'positive' refers to state where covariate has a zero effect.
  # . 'negative' refers to state where covariate has a non-zero effect.
  # Thus, true-positive refers to state where covariate has a zero effect
  # and the selection methods excludes covariate.
  # True-negative refers to state where covariate has a non-zero effect
  # and the selection method does not exclude covariate.
  # False-positive : covariate has a zero effect and is not excluded.
  # False-negative : covariate has non-zero effect and is not excluded.

  # True positive : an un-associated variable was excluded
  tp_selection_count_G <- sum((covariate_identifier_G[1,] == covariate_identifier_G[2,]) & (covariate_identifier_G[1,] == 0))
  # True negative : an associated variable was not excluded
  tn_selection_count_G <- sum((covariate_identifier_G[1,] == covariate_identifier_G[2,]) & (covariate_identifier_G[1,] != 0))
  # False positive : an associated variable was excluded
  fp_selection_count_G <- sum((covariate_identifier_G[1,] != covariate_identifier_G[2,]) & (covariate_identifier_G[1,] != 0))
  # False negative : an un-associated variable was not excluded
  fn_selection_count_G <- sum((covariate_identifier_G[1,] != covariate_identifier_G[2,]) & (covariate_identifier_G[1,] == 0))


  #TODO
  tp_selection_rate_G <- tp_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  tn_selection_rate_G <- tn_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  fp_selection_rate_G <- fp_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)
  fn_selection_rate_G <- fn_selection_count_G #/ n_G_attr #/(nonzero_controls+1e-4)


  selection_confounder_identifier <- NaN
  if (dgp == "houseprices") {
    # Compute share of individual confounder that was excluded from model
    selection_confounder_identifier <- selection_identifier
    selection_confounder_identifier[is.na(selection_confounder_identifier)] <- 0
    selection_confounder_identifier <- selection_confounder_identifier[1]
    #selection_confounder_identifier <- covariate_identifier_G[1]
  }

  # Compute metrics of interest (see below) with test data.
  y_test <- data_test$y
  D_test <- data_test$D

  # Assemble different regressor matrices
  Z_test <- data.matrix(data_test[colnames_covariates])

  Z_test_selected <- Z_test[, selection_identifier[!is.na(selection_identifier)]]
  X_test_selected <- data.matrix(cbind(D_test, Z_test_selected))

  # Estimate treatment effect by OLS provided with selected covariates
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
              selection_confounder_identifier=selection_confounder_identifier))
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
