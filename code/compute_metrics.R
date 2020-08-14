### CompStat project --- Max Schäfer
# This script contains the function computing all metrics used
# to evaluate the variable selection methods


compute_metrics <- function(n=200, n_F_attr=30, n_G_attr=30, n_H_attr=30, treatment_effect=.5,
                                    beta_GD_size=.5, beta_GY_size=.25, beta_H_size=1,
                                    beta_F_size=1,
                                    unconfoundedness_rate=.9, selection_method="simple") {


  # Draw data once and split into train and test (*)
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
    #TODO
    out <- simple_select_covariates(data_train=data_train)
    selection_identifier <- out$selected_covars_one
    select_treatment_identifier <- out$select_treatment_effect
  }
  if (selection_method=="double") {
    selection_identifier <- double_select_covariates(data_train=data_train)
    select_treatment_identifier <- NaN
  }


  if (length(true_covariate_identifier) != length(selection_identifier)) {
    print("Warning: check how true_covariate_identifier vector is initialized in data_generating_process.R")
  }
  # Collect information which potential confounders were selected
  covariate_identifier <- data.frame(rbind(true_covariate_identifier,
                                           selection_identifier))

  #names(covariate_identifier) <- colnames_GH
  names(covariate_identifier) <- colname_G
  covariate_identifier_G <- covariate_identifier[colname_G]
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
  X_test <- data.matrix(cbind(D_test, G_test, H_test))
  Z_test <- data.matrix(cbind(G_test, H_test))

  X_test <- data.matrix(cbind(D_test, G_test))
  Z_test <- data.matrix(G_test)

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
