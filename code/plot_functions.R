### CompStat project --- Max Schäfer


lasso_coef_shrink_plot <- function(df, lambda_min=NaN, figure_name="", title="") {
  # Plot all lasso coefficients for a grid of lambda values
  # Distinguish between zero and non-zero effect controls

  cols <- names(df)[1:(length(names(df))-1)]

  # Initialize plot.
  plot <- ggplot(data=df, aes(x=lambda_grid))

  # Loop through each beta and plot its value with increasing penalty.
  for (name in cols) {
    random_transparency <- runif(n=1, min=.5, max=1) #

    #if (startsWith(name, "F")) {
    #  plot <- plot +
    #          geom_line(aes_string(y=name), size=1, col="#a89932", alpha=random_transparency)
    #}
    #else if (startsWith(name, "G")) {
    #  plot <- plot +
    #          geom_line(aes_string(y=name), size=1, col="#f74376", alpha=random_transparency)
    #}
    if (startsWith(name, "G") & endsWith(name, "p")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="green", alpha=.5*random_transparency)
    }
    else if (startsWith(name, "G") & endsWith(name, "n")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="red", alpha=random_transparency)
    }
    #else if (startsWith(name, "H")) {
    #  plot <- plot +
    #          geom_line(aes_string(y=name), size=1, col="#6143f7", alpha=random_transparency)
    #}
    else if (startsWith(name, "D")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#586e5e")
    }
  }
  plot <- plot +
          geom_line(aes(y=0), lty="dashed") +
          ggtitle(figure_name, subtitle=title) +
          ylab("Coefficient Value") +
          xlab("Lambda") +
          theme_grey(base_size=15)

  if (is.na(lambda_min) == FALSE) {
    plot <- plot + geom_vline(xintercept=lambda_min, linetype="dashed",
                              color="black", alpha=.8, size=1)
  }

  return(plot)
}






# Variable selection plot (TP, TN, FP, FN)


variable_selection_rate_plot <- function(df, selection_method, subtitle="", xlab="", ylim=NaN) {

  plot <- ggplot(data=df, aes(x=x_axis))
  plot <- plot +
          geom_line(aes_string(y="tp_selection_rate_G"), size=1, col="blue", alpha=.9) +
          geom_line(aes_string(y="tn_selection_rate_G"), size=1, col="blue", alpha=.5) +
          geom_line(aes_string(y="fp_selection_rate_G"), size=1, col="red", alpha=.9) +
          geom_line(aes_string(y="fn_selection_rate_G"), size=1, col="red", alpha=.5) +
          geom_line(aes(y=0), lty="dashed") +
          labs(title="TP, TN, FP, FN rates",
                     subtitle=subtitle,
                     caption="Legend: Dark blue shows True-Positive, Light blue is True-Negative, Dark red is False-Positive, Light red is False-Negative",
                     x=xlab, y="") +
          theme_grey(base_size=15)

  if (is.na(ylim[1]) == FALSE) {
    plot <- plot + ylim(ylim)
  }

  return(plot)
}


confounding_bias_plot <- function(df, selection_method, subtitle="", xlab="", ylim=NaN) {

  plot <- ggplot(data=df, aes(x=x_axis))
  plot <- plot +
          geom_line(aes_string(y="root_mean_squared_bias"), size=1, col="#f5ce42", alpha=1) +
          geom_line(aes_string(y="mean_abs_dev_bias"), size=1, col="red", alpha=.7) +
          geom_line(aes_string(y="median_bias"), size=1, col="#bf42f5", alpha=1) +
          geom_line(aes(y=0), lty="dashed") +
          labs(title="Bias of Treatment Coefficient",
                     subtitle=subtitle,
                     caption="Legend: Purple shows Median Bias, Orange is Root Mean Squared Bias, Red is Average Absolute Deviation",
                     x=xlab, y="") +
          theme_grey(base_size=15)

  if (is.na(ylim[1]) == FALSE) {
    plot <- plot + ylim(ylim)
  }

  return(plot)
}



purchasing_power_exclusion_plot <- function(df, selection_method, subtitle="", xlab="", ylim=NaN) {

  plot <- ggplot(data=df, aes(x=x_axis))
  plot <- plot +
          geom_line(aes_string(y="selection_confounder_identifier"), size=1, col="red", alpha=1) +
          geom_line(aes(y=0), lty="dashed") +
          labs(title="Selection of the Confounder 'Purchasing Power'",
                     subtitle=subtitle,
                     caption="Source: ",
                     x=xlab, y="")

  if (is.na(ylim[1]) == FALSE) {
    plot <- plot + ylim(ylim)
  }

  return(plot)
}




sim_plot_wrapper <- function(dgp, R, iter_over, sim_parameter_vec,
                             selection_method,
                             n, n_F_attr,
                             n_G_attr, n_H_attr,
                             treatment_effect, beta_GD_size,
                             beta_GY_size, beta_F_size, beta_H_size,
                             nonzero_controls,
                             ylim_select=NaN, ylim_bias=NaN) {


  simulation_results <- simulation_wrapper(dgp=dgp, R=R, iter_over=iter_over,
                              sim_parameter_vec=sim_parameter_vec,
                              selection_method=selection_method,
                              n=n, n_F_attr=n_F_attr,
                              n_G_attr=n_G_attr, n_H_attr=n_H_attr,
                              treatment_effect=treatment_effect,
                              beta_GD_size=beta_GD_size,
                              beta_GY_size=beta_GY_size,
                              beta_F_size=beta_F_size,
                              beta_H_size=beta_H_size,
                              nonzero_controls=nonzero_controls)

  root_mean_squared_bias <- simulation_results$root_mean_squared_bias_vec
  mean_abs_dev_bias <- simulation_results$mean_abs_dev_bias_vec
  median_bias <- simulation_results$median_bias_vec
  tp_selection_rate_G <- simulation_results$tp_selection_rate_G_vec
  tn_selection_rate_G <- simulation_results$tn_selection_rate_G_vec
  fp_selection_rate_G <- simulation_results$fp_selection_rate_G_vec
  fn_selection_rate_G <- simulation_results$fn_selection_rate_G_vec
  selection_confounder_identifier <- simulation_results$selection_confounder_identifier_vec

  # Visualization
  x_axis <- sim_parameter_vec
  # Metrics --- variable selection
  df <- data.frame(cbind(x_axis,
                         tp_selection_rate_G, tn_selection_rate_G,
                         fp_selection_rate_G, fn_selection_rate_G))
  # Plot
  p_variable_select <- variable_selection_rate_plot(df=df, selection_method=selection_method, subtitle=str_c(selection_method, "-selection"), xlab=iter_over, ylim=ylim_select)

  # Show exclusion rate of particular confounder, default: 'purchasing power'
  df <- data.frame(cbind(x_axis, selection_confounder_identifier))
  p_purchasing_power_exclusion <- purchasing_power_exclusion_plot(df=df, subtitle=str_c(selection_method, "-selection"), xlab=iter_over, ylim=c(-.1, 1.1))

  # Metrics --- confounder bias
  df <- data.frame(cbind(x_axis, root_mean_squared_bias, mean_abs_dev_bias, median_bias))
  # Plot
  p_confounder_bias <- confounding_bias_plot(df=df, selection_method=selection_method, subtitle=str_c(selection_method, "-selection"), xlab=iter_over, ylim=ylim_bias)

  return(list(p_variable_select=p_variable_select,
              p_confounder_bias=p_confounder_bias,
              p_purchasing_power_exclusion=p_purchasing_power_exclusion))
}





produce_lasso_coef_shrink_plot <- function(dgp, n, n_F_attr=NaN, n_G_attr, n_H_attr=NaN,
                                           treatment_effect, beta_GD_size, beta_GY_size,
                                           beta_F_size=.5, beta_H_size=.5,
                                           nonzero_controls,
                                           lambda_grid=seq(0, 1, by=.02),
                                           title="") {


  # Update column names for identifying regressors
  colnames_confounders <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))

  if (dgp == "A") {
    colnames_covariates <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
    colnumbers_covariates <- seq(from=1, to=n_G_attr, by=1)
    # Draw dataset which is then split into training and test
    data_set <- generate_data_A(n=n,
                                n_G_attr=n_G_attr,
                                treatment_effect=treatment_effect,
                                beta_GD_size=beta_GD_size,
                                beta_GY_size=beta_GY_size,
                                nonzero_controls=nonzero_controls)
  }
  if (dgp == "B") {
    colnames_covariates <- c(str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1)),
                             str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1)),
                             str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1)))

    colnumbers_covariates <- seq(from=1, to=(n_F_attr + n_G_attr + n_H_attr), by=1)
    # Draw dataset which is then split into training and test
    data_set <- generate_data_B(n=n,
                                n_F_attr=n_F_attr,
                                n_G_attr=n_G_attr,
                                n_H_attr=n_H_attr,
                                treatment_effect=treatment_effect,
                                beta_GD_size=beta_GD_size,
                                beta_GY_size=beta_GY_size,
                                beta_H_size=beta_H_size,
                                beta_F_size=beta_F_size,
                                nonzero_controls=nonzero_controls)
  }
  if (dgp == "houseprices") {
    colnames_covariates <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
    colnumbers_covariates <- seq(from=1, to=n_G_attr, by=1)
    # Draw dataset which is then split into training and test
    data_set <- generate_data_houseprices_dgp0(n=n, n_G_attr=n_G_attr, corr_G=0,
                                               treatment_effect=-7, beta_GY_inflator=1)
  }



  data <- data_set$data
  # Shuffle data before train-test-split to ensure no unintended structure exists
  # due to the order of observations
  random_rows <- sample(nrow(data))
  data <- data[random_rows, ]

  # Split data into training and test sets of equal size
  data_train <- data[1:(dim(data)[1]/2), ]
  data_test <- data[(dim(data)[1]/2+1):dim(data)[1], ]

  y_train <- data_train$y
  D_train <- data_train$D
  # Assemble different regressor matrices
  Z_train <- data.matrix(data_train[colnames_covariates])
  X_train <- data.matrix(cbind(D_train, Z_train))



  # Retrieve true covariate identifier
  true_covariate_identifier <- data_set$true_covariate_identifier
  names(true_covariate_identifier) <- colnames_covariates
  true_covariate_identifier[true_covariate_identifier > 0] <- "n"
  true_covariate_identifier[true_covariate_identifier != "n"] <- "p"

  # Change column name to identify sparse covariates
  colnames_two_effect <- str_c(colnames_covariates, true_covariate_identifier)
  colnames_one_effect <- c("D", colnames_two_effect)



  # Variable selection in "first" and "second" stage

  # Grid of lasso penalties
  lambda_grid <- lambda_grid

  # 1. Use Lasso to shrink/select covariates given their association with outcome Y
  lasso_one <- glmnet(X_train, y_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

  # 2. Use Lasso to shrink/select covariates given their association with tratment D
  lasso_two <- glmnet(Z_train, D_train, alpha=1, lambda=lambda_grid, intercept=FALSE)

  # Extract estimates from both models
  beta_hats_lasso_one <- t(matrix(coef(lasso_one, s=lambda_grid)[-1,], nrow=dim(X_train)[2], ncol=length(lambda_grid)))
  beta_hats_lasso_two <- t(matrix(coef(lasso_two, s=lambda_grid)[-1,], nrow=dim(Z_train)[2], ncol=length(lambda_grid)))


  # MSE-optimal penalty term for lasso
  lasso_one_cv <- cv.glmnet(X_train, y_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lambda_min_one <- lasso_one_cv$lambda.min

  lasso_two_cv <- cv.glmnet(Z_train, D_train, alpha=1, intercept=FALSE,
                            type.measure="mse", nfolds=10)
  lambda_min_two <- lasso_two_cv$lambda.min

  # 1.
  df <- data.frame(cbind(beta_hats_lasso_one, lambda_grid))
  names(df) <- colnames_one_effect
  p1 <- lasso_coef_shrink_plot(df=df, lambda_min=lambda_min_one, figure_name="", title=str_c("LASSO regression of y on X", title))


  # 2.
  df <- data.frame(cbind(beta_hats_lasso_two, lambda_grid))
  names(df) <- colnames_two_effect
  p2 <- lasso_coef_shrink_plot(df=df, lambda_min=lambda_min_two, title=str_c("LASSO regression of D on X", title))

  return(list(p1=p1, p2=p2))
}




produce_bias_hist <- function(dgp, R, n, n_F_attr, n_G_attr, n_H_attr,
                              treatment_effect, beta_GD_size, beta_GY_size,
                              beta_F_size, beta_H_size, nonzero_controls,
                              binwidth=.04, title="",
                              xlim=c(-.4, .4), ylim=c(0, 10)) {

  selection_method <- "simple"
  simulation_results <- simulation_wrapper(dgp=dgp, R=R, iter_over=NaN,
                                           sim_parameter_vec=rep(NA, 3),
                                           selection_method=selection_method,
                                           n=n,
                                           n_F_attr=n_F_attr, n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           treatment_effect=treatment_effect,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_F_size=beta_F_size,
                                           beta_H_size=beta_H_size,
                                           nonzero_controls=nonzero_controls)

  confounder_bias_vec_simple <- simulation_results$confounder_bias_vec
  tp_selection_rate_G_simple <- simulation_results$tp_selection_rate_G_vec
  fp_selection_rate_G_simple <- simulation_results$fp_selection_rate_G_vec


  selection_method <- "double"
  simulation_results <- simulation_wrapper(dgp=dgp, R=R, iter_over=NaN,
                                           sim_parameter_vec=rep(NA, 3),
                                           selection_method=selection_method,
                                           n=n,
                                           n_F_attr=n_F_attr, n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           treatment_effect=treatment_effect,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_F_size=beta_F_size,
                                           beta_H_size=beta_H_size,
                                           nonzero_controls=nonzero_controls)

  confounder_bias_vec_double <- simulation_results$confounder_bias_vec
  tp_selection_rate_G_double <- simulation_results$tp_selection_rate_G_vec
  fp_selection_rate_G_double <- simulation_results$fp_selection_rate_G_vec


  df <- data.frame(confounder_bias_vec_simple)
  hist_simple <- ggplot(data=df, aes_string(x="confounder_bias_vec_simple")) +
                  geom_histogram(aes(y=..density..), binwidth=binwidth, colour="blue", fill="white") +
                  geom_density(alpha=.2, fill="#FF6666") +
                  geom_vline(xintercept=0, linetype="dashed", color="black", alpha=.8, size=1) +
                  labs(title=str_c("Bias of treatment effect estimate", title),
                       subtitle="simple-LASSO",
                       caption=str_c("Notes: Average False-Positive count is ", round(mean(fp_selection_rate_G_simple), 2)),
                       x="Confounder bias", y="Density") +
                  xlim(xlim) +
                  ylim(ylim) +
                  theme_grey(base_size=15)


  df <- data.frame(confounder_bias_vec_double)
  hist_double <- ggplot(data=df, aes_string(x="confounder_bias_vec_double")) +
                  geom_histogram(aes(y=..density..), binwidth=binwidth, colour="blue", fill="white") +
                  geom_density(alpha=.2, fill="#FF6666") +
                  geom_vline(xintercept=0, linetype="dashed", color="black", alpha=.8, size=1) +
                  labs(title=str_c("Bias of treatment effect estimate", title),
                       subtitle="double-LASSO",
                       caption=str_c("Notes: Average False-Positive count is ", round(mean(fp_selection_rate_G_double), 2)),
                       x="Confounder bias", y="Density") +
                  xlim(xlim) +
                  ylim(ylim) +
                  theme_grey(base_size=15)

  return(list(hist_simple=hist_simple, hist_double=hist_double))
}





































#
