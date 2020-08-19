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
          xlab("Lambda")

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
                     caption="Source: ",
                     x=xlab, y="")

  #geom_line(aes_string(y=str_c("tp_selection_rate_G_", selection_method)), size=1, col="blue", alpha=.9) +
  #geom_line(aes_string(y=str_c("tn_selection_rate_G_", selection_method)), size=1, col="blue", alpha=.5) +
  #geom_line(aes_string(y=str_c("fp_selection_rate_G_", selection_method)), size=1, col="red", alpha=.9) +
  #geom_line(aes_string(y=str_c("fn_selection_rate_G_", selection_method)), size=1, col="red", alpha=.5) +

  if (is.na(ylim) == FALSE) {
    plot <- plot + ylim(ylim)
  }

  return(plot)
}


confounding_bias_plot <- function(df, selection_method, subtitle="", xlab="", ylim=NaN) {

  plot <- ggplot(data=df, aes(x=x_axis))
  plot <- plot +
          geom_line(aes_string(y="mean_squared_bias"), size=1, col="red", alpha=1) +
          geom_line(aes_string(y="mean_abs_dev_bias"), size=1, col="#575757", alpha=.6) +
          geom_line(aes_string(y="median_bias"), size=1, col="#575757", alpha=1) +
          geom_line(aes(y=0), lty="dashed") +
          labs(title="Bias of Treatment Coefficient",
                     subtitle=subtitle,
                     caption="Source: ",
                     x=xlab, y="")

  if (is.na(ylim) == FALSE) {
    plot <- plot + ylim(ylim)
  }
  #geom_line(aes_string(y=str_c("mean_squared_bias_", selection_method)), size=1, col="red", alpha=1) +
  #geom_line(aes_string(y=str_c("mean_abs_dev_bias_", selection_method)), size=1, col="#575757", alpha=.6) +
  #geom_line(aes_string(y=str_c("median_bias_", selection_method)), size=1, col="#575757", alpha=1) +

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

  if (is.na(ylim) == FALSE) {
    plot <- plot + ylim(ylim)
  }
  #geom_line(aes_string(y=str_c("mean_squared_bias_", selection_method)), size=1, col="red", alpha=1) +
  #geom_line(aes_string(y=str_c("mean_abs_dev_bias_", selection_method)), size=1, col="#575757", alpha=.6) +
  #geom_line(aes_string(y=str_c("median_bias_", selection_method)), size=1, col="#575757", alpha=1) +

  return(plot)
}




sim_plot_wrapper <- function(dgp, R, iter_over, sim_parameter_vec,
                             selection_method,
                             n, n_F_attr,
                             n_G_attr, n_H_attr,
                             treatment_effect, beta_GD_size,
                             beta_GY_size, beta_F_size, nonzero_controls,
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
                              nonzero_controls=nonzero_controls)

  mean_squared_bias <- simulation_results$mean_squared_bias_vec
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
  df <- data.frame(cbind(x_axis, mean_squared_bias, mean_abs_dev_bias, median_bias))
  # Plot
  p_confounder_bias <- confounding_bias_plot(df=df, selection_method=selection_method, subtitle=str_c(selection_method, "-selection"), xlab=iter_over, ylim=ylim_bias)

  return(list(p_variable_select=p_variable_select,
              p_confounder_bias=p_confounder_bias,
              p_purchasing_power_exclusion=p_purchasing_power_exclusion))
}















#
