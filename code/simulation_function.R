### CompStat project --- Max Schäfer


# Wrapper for simulation function
simulation_wrapper <- function(dgp, R, iter_over, sim_parameter_vec,
                               selection_method,
                               n, n_F_attr, n_G_attr, n_H_attr, corr_G=0,
                               treatment_effect, beta_GD_size,
                               beta_GY_size, beta_F_size, beta_H_size,
                               nonzero_controls) {

  # Housekeeping
  colnames_confounders <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
  colnames_covariates <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))

  if (dgp == "B") {
    colnames_covariates <- c(str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1)),
                             str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1)),
                             str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1)))
  }
  if (dgp == "houseprices") {
    #colnames_covariates <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
    #TODO above for dgp-a-houseprices, below for houseprices (in line with -b)
    colnames_covariates <- c(str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1)),
                             str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1)),
                             str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1)))
  }


  # Simulation with R runs WITHOUT varying over a parameter.
  if (is.na(iter_over) == TRUE) {

    # Initialize containers
    confounder_bias_vec <- rep(NA, R)
    test_MSE_vec <- rep(NA, R)
    tp_selection_rate_G_vec <- rep(NA, R)
    tn_selection_rate_G_vec <- rep(NA, R)
    fp_selection_rate_G_vec <- rep(NA, R)
    fn_selection_rate_G_vec <- rep(NA, R)
    selection_confounder_identifier_vec <- rep(NA, R)

    for (i in 1:R) {
      if (dgp == "A") {
        data_set <- generate_data_A(n=n,
                                    n_G_attr=n_G_attr,
                                    corr_G=0,
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
                                    corr_G=0,
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
      true_covariate_identifier <- data_set$true_covariate_identifier

      simulation_results <- compute_metrics(data=data,
                                            true_covariate_identifier=true_covariate_identifier,
                                            selection_method=selection_method,
                                            colnames_covariates=colnames_covariates,
                                            colnames_confounders=colnames_confounders,
                                            treatment_effect=treatment_effect)

      # Metrics --- confounder bias
      confounder_bias_vec[i] <- simulation_results$confounder_bias
      # Metrics --- prediction
      test_MSE_vec[i] <- simulation_results$test_MSE
      # Metrics --- variable selection
      tp_selection_rate_G_vec[i] <- simulation_results$tp_selection_rate_G
      tn_selection_rate_G_vec[i] <- simulation_results$tn_selection_rate_G
      fp_selection_rate_G_vec[i] <- simulation_results$fp_selection_rate_G
      fn_selection_rate_G_vec[i] <- simulation_results$fn_selection_rate_G

      selection_confounder_identifier_vec[i] <- simulation_results$selection_confounder_identifier
    }

    return(list(confounder_bias_vec=confounder_bias_vec,
                test_MSE_vec=test_MSE_vec,
                tp_selection_rate_G_vec=tp_selection_rate_G_vec,
                tn_selection_rate_G_vec=tn_selection_rate_G_vec,
                fp_selection_rate_G_vec=fp_selection_rate_G_vec,
                fn_selection_rate_G_vec=fn_selection_rate_G_vec,
                selection_confounder_identifier_vec=selection_confounder_identifier_vec))
  }



  # Simulation with R runs WHILE varying over a parameter.
  if (is.na(iter_over) == FALSE) {

    # Initialize containers
    root_mean_squared_bias_vec <- rep(NA, length(sim_parameter_vec))
    mean_abs_dev_bias_vec <- rep(NA, length(sim_parameter_vec))
    median_bias_vec <- rep(NA, length(sim_parameter_vec))
    mean_test_MSE_vec <- rep(NA, length(sim_parameter_vec))
    tp_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))
    tn_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))
    fp_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))
    fn_selection_rate_G_vec <- rep(NA, length(sim_parameter_vec))

    selection_confounder_identifier_vec <- rep(NA, length(sim_parameter_vec))



    i <- 1
    start_time <- Sys.time()

    for (SIM_PARAMETER in sim_parameter_vec) {

      if (iter_over=="nonzero_controls") {
        # Draw data once, and then split into train and test (*)
        if (dgp == "A") {
          data_set <- generate_data_A(n=n,
                                      n_G_attr=n_G_attr,
                                      corr_G=corr_G,
                                      treatment_effect=treatment_effect,
                                      beta_GD_size=beta_GD_size,
                                      beta_GY_size=beta_GY_size,
                                      nonzero_controls=SIM_PARAMETER)
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
                                      nonzero_controls=SIM_PARAMETER)
        }

        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }

      if (iter_over=="beta_GY_size") {
        if (dgp == "A") {
          data_set <- generate_data_A(n=n,
                                      n_G_attr=n_G_attr,
                                      corr_G=corr_G,
                                      treatment_effect=treatment_effect,
                                      beta_GD_size=beta_GD_size,
                                      beta_GY_size=SIM_PARAMETER,
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
                                      beta_GY_size=SIM_PARAMETER,
                                      beta_H_size=beta_H_size,
                                      beta_F_size=beta_F_size,
                                      nonzero_controls=nonzero_controls)
        }
        if (dgp == "houseprices") {
          #data_set <- generate_data_houseprices_dgp0(n=n, n_G_attr=n_G_attr,
          #                                           corr_G=corr_G,
          #                                           beta_GY_inflator=SIM_PARAMETER)
          data_set <- generate_data_houseprices(n=n, n_G_attr=n_G_attr,
                                                n_F_attr=n_F_attr, n_H_attr=n_H_attr,
                                                corr_G=corr_G,
                                                beta_GY_inflator=SIM_PARAMETER)
        }

        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }
      if (iter_over=="beta_GD_size") {
        if (dgp == "A") {
          data_set <- generate_data_A(n=n,
                                      n_G_attr=n_G_attr,
                                      corr_G=corr_G,
                                      treatment_effect=treatment_effect,
                                      beta_GD_size=SIM_PARAMETER,
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
                                      beta_GD_size=SIM_PARAMETER,
                                      beta_GY_size=beta_GY_size,
                                      beta_H_size=beta_H_size,
                                      beta_F_size=beta_F_size,
                                      nonzero_controls=nonzero_controls)
        }
        if (dgp == "houseprices") {
          data_set <- generate_data_houseprices_dgp0(n=n, n_G_attr=n_G_attr,
                                                     corr_G=corr_G,
                                                     beta_GY_inflator=beta_GY_inflator)
        }
        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }
      if (iter_over=="n_G_attr") {
        colnames_confounders <- str_c(rep("G_", SIM_PARAMETER), seq(from=1, to=SIM_PARAMETER, by=1))
        colnames_covariates <- str_c(rep("G_", SIM_PARAMETER), seq(from=1, to=SIM_PARAMETER, by=1))

        if (dgp == "A") {
          data_set <- generate_data_A(n=n,
                                      n_G_attr=SIM_PARAMETER,
                                      corr_G=corr_G,
                                      treatment_effect=treatment_effect,
                                      beta_GD_size=beta_GD_size,
                                      beta_GY_size=beta_GY_size,
                                      nonzero_controls=nonzero_controls)
        }
        if (dgp == "B") {
          colnames_covariates <- c(str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1)),
                                   str_c(rep("G_", SIM_PARAMETER), seq(from=1, to=SIM_PARAMETER, by=1)),
                                   str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1)))
          data_set <- generate_data_B(n=n,
                                      n_G_attr=SIM_PARAMETER,
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
          data_set <- generate_data_houseprices_dgp0(n=n, n_G_attr=SIM_PARAMETER,
                                                     corr_G=corr_G,
                                                     beta_GY_inflator=beta_GY_size)
        }

        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }


      ###


      if (iter_over=="n_F_attr") {
        colname_F <- str_c(rep("F_", SIM_PARAMETER), seq(from=1, to=SIM_PARAMETER, by=1))
        if (dgp == "B") {
          data_set <- generate_data_B(n=n,
                                      n_G_attr=n_G_attr,
                                      n_F_attr=SIM_PARAMETER,
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
        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }
      if (iter_over=="beta_F_size") {
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
                                      beta_F_size=SIM_PARAMETER,
                                      nonzero_controls=nonzero_controls)
        }
        if (dgp == "houseprices") {
          data_set <- generate_data_houseprices_dgp0(n=n, n_G_attr=n_G_attr,
                                                     corr_G=corr_G,
                                                     beta_GY_inflator=beta_GY_size)
        }
        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }
      if (iter_over=="corr_confounders") {
        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }
      if (iter_over=="treatment_effect") {
        sim_results_vec <- replicate(R, compute_metrics(data=data_set$data,
                                              true_covariate_identifier=data_set$true_covariate_identifier,
                                              selection_method=selection_method,
                                              colnames_covariates=colnames_covariates,
                                              colnames_confounders=colnames_confounders,
                                              treatment_effect=treatment_effect))
      }

      # Metrics --- confounder bias
      root_mean_squared_bias_vec[i] <- sqrt(mean((unlist(sim_results_vec[1, 1:R]))**2))
      mean_abs_dev_bias_vec[i] <- mean(abs(unlist(sim_results_vec[1, 1:R])))
      median_bias_vec[i] <- median(unlist(sim_results_vec[1, 1:R]))
      # Metrics --- prediction
      mean_test_MSE_vec[i] <- mean(unlist(sim_results_vec[2, 1:R]))
      # Metrics --- variable selection
      tp_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[3, 1:R]))
      tn_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[4, 1:R]))
      fp_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[5, 1:R]))
      fn_selection_rate_G_vec[i] <- mean(unlist(sim_results_vec[6, 1:R]))

      selection_confounder_identifier_vec[i] <- mean(unlist(sim_results_vec[7, 1:R]))

      end_time <- Sys.time()
      if (i == 1) {
        delta <- end_time - start_time
        #cat(" ... this takes approximately", round(R*delta, 2), "seconds.")
      }
      i <- i+1
    }
    return(list(root_mean_squared_bias_vec=root_mean_squared_bias_vec,
                mean_abs_dev_bias_vec=mean_abs_dev_bias_vec,
                median_bias_vec=median_bias_vec,
                mean_test_MSE_vec=mean_test_MSE_vec,
                tp_selection_rate_G_vec=tp_selection_rate_G_vec,
                tn_selection_rate_G_vec=tn_selection_rate_G_vec,
                fp_selection_rate_G_vec=fp_selection_rate_G_vec,
                fn_selection_rate_G_vec=fn_selection_rate_G_vec,
                selection_confounder_identifier_vec=selection_confounder_identifier_vec))
  }
}
