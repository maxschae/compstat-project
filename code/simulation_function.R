### CompStat project --- Max Schäfer

# Wrapper for simulation function
simulation_wrapper <- function(R, iter_over, sim_parameter_vec,
                               selection_method,
                               n, n_F_attr, n_G_attr, n_H_attr,
                               treatment_effect, beta_GD_size,
                               beta_GY_size, beta_F_size,
                               unconfoundedness_rate) {

  #
  cat("SIMULATION STARTED: vary over", iter_over, "with", R,
      "repititions for each value")

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
  start_time <- Sys.time()

  for (SIM_PARAMETER in sim_parameter_vec) {

    if (iter_over=="unconfoundedness_rate") {
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=n_F_attr,
                                           n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_H_size=beta_H_size,
                                           beta_F_size=beta_F_size,
                                           treatment_effect=treatment_effect,
                                           unconfoundedness_rate=SIM_PARAMETER,
                                           selection_method=selection_method))
    }
    if (iter_over=="treatment_effect") {
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=n_F_attr,
                                           n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_H_size=beta_H_size,
                                           beta_F_size=beta_F_size,
                                           treatment_effect=SIM_PARAMETER,
                                           unconfoundedness_rate=unconfoundedness_rate,
                                           selection_method=selection_method))
    }
    if (iter_over=="beta_GY_size") {
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=n_F_attr,
                                           n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=SIM_PARAMETER,
                                           beta_H_size=beta_H_size,
                                           beta_F_size=beta_F_size,
                                           treatment_effect=treatment_effect,
                                           unconfoundedness_rate=unconfoundedness_rate,
                                           selection_method=selection_method))
    }
    if (iter_over=="beta_GD_size") {
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=n_F_attr,
                                           n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=SIM_PARAMETER,
                                           beta_GY_size=beta_GY_size,
                                           beta_H_size=beta_H_size,
                                           beta_F_size=beta_F_size,
                                           treatment_effect=treatment_effect,
                                           unconfoundedness_rate=unconfoundedness_rate,
                                           selection_method=selection_method))
    }
    if (iter_over=="beta_H_size") {
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=n_F_attr,
                                           n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_H_size=SIM_PARAMETER,
                                           beta_F_size=beta_F_size,
                                           treatment_effect=treatment_effect,
                                           unconfoundedness_rate=unconfoundedness_rate,
                                           selection_method=selection_method))
    }
    if (iter_over=="n_G_attr") {
      colname_G <- str_c(rep("G_", SIM_PARAMETER), seq(from=1, to=SIM_PARAMETER, by=1))
      colnames_GH <- c(colname_G, colname_H)
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=n_F_attr,
                                           n_G_attr=SIM_PARAMETER,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_H_size=beta_H_size,
                                           beta_F_size=beta_F_size,
                                           treatment_effect=treatment_effect,
                                           unconfoundedness_rate=unconfoundedness_rate,
                                           selection_method=selection_method))
    }
    if (iter_over=="n_F_attr") {
      colname_F <- str_c(rep("F_", SIM_PARAMETER), seq(from=1, to=SIM_PARAMETER, by=1))
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=SIM_PARAMETER,
                                           n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_H_size=beta_H_size,
                                           beta_F_size=beta_F_size,
                                           treatment_effect=treatment_effect,
                                           unconfoundedness_rate=unconfoundedness_rate,
                                           selection_method=selection_method))
    }
    if (iter_over=="beta_F_size") {
      sim_results_vec <- replicate(R, compute_metrics(n=n,
                                           n_F_attr=n_F_attr,
                                           n_G_attr=n_G_attr,
                                           n_H_attr=n_H_attr,
                                           beta_GD_size=beta_GD_size,
                                           beta_GY_size=beta_GY_size,
                                           beta_H_size=beta_H_size,
                                           beta_F_size=SIM_PARAMETER,
                                           treatment_effect=treatment_effect,
                                           unconfoundedness_rate=unconfoundedness_rate,
                                           selection_method=selection_method))
    }


    # Metrics --- confounder bias
    #confounder_bias_mat[i, ] <- unlist(sim_results_vec[1, 1:R]) #TODO
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

    end_time <- Sys.time()
    if (i == 1) {
      delta <- end_time - start_time
      cat(" ... this takes approximately", round(R*delta, 2), "seconds.")
    }
    i <- i+1
  }


  return(list(mean_squared_bias_vec=mean_squared_bias_vec,
              mean_abs_dev_bias_vec=mean_abs_dev_bias_vec,
              median_bias_vec=median_bias_vec,
              mean_test_MSE_vec=mean_test_MSE_vec,
              tp_selection_rate_G_vec=tp_selection_rate_G_vec,
              tn_selection_rate_G_vec=tn_selection_rate_G_vec,
              fp_selection_rate_G_vec=fp_selection_rate_G_vec,
              fn_selection_rate_G_vec=fn_selection_rate_G_vec,
              select_treatment_identifier_vec=select_treatment_identifier_vec))
}
