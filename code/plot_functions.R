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

    if (startsWith(name, "F")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#a89932", alpha=random_transparency)
    }
    #else if (startsWith(name, "G")) {
    #  plot <- plot +
    #          geom_line(aes_string(y=name), size=1, col="#f74376", alpha=random_transparency)
    #}
    else if (startsWith(name, "G") & endsWith(name, "p")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="green", alpha=.5*random_transparency)
    }
    else if (startsWith(name, "G") & endsWith(name, "n")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="red", alpha=random_transparency)
    }
    else if (startsWith(name, "H")) {
      plot <- plot +
              geom_line(aes_string(y=name), size=1, col="#6143f7", alpha=random_transparency)
    }
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


variable_selection_rate_plot <- function(df, subtitle="", xlab="", ylim=NaN) {

  plot <- ggplot(data=df, aes(x=x_axis))
  plot <- plot +
          geom_line(aes_string(y="tp_selection_rate_G_simple"), size=1, col="blue", alpha=.9) +
          geom_line(aes_string(y="tn_selection_rate_G_simple"), size=1, col="blue", alpha=.5) +
          geom_line(aes_string(y="fp_selection_rate_G_simple"), size=1, col="red", alpha=.9) +
          geom_line(aes_string(y="fn_selection_rate_G_simple"), size=1, col="red", alpha=.5) +
          geom_line(aes(y=0), lty="dashed") +
          ggtitle("TP, TN, FP, FN rates", subtitle=subtitle) +
          ylab("") +
          xlab(xlab)

  if (is.na(ylim) == FALSE) {
    plot <- plot + ylim(ylim)
  }

  return(plot)
}


confounding_bias_plot <- function(df, subtitle="", xlab="", ylim=c(-.1, .3)) {

  plot <- ggplot(data=df, aes(x=x_axis))
  plot <- plot +
          geom_line(aes_string(y="mean_squared_bias_simple"), size=1, col="red", alpha=1) +
          geom_line(aes_string(y="mean_abs_dev_bias_simple"), size=1, col="#575757", alpha=.6) +
          geom_line(aes_string(y="median_bias_simple"), size=1, col="#575757", alpha=1) +
          geom_line(aes(y=0), lty="dashed") +
          ylim(ylim) +
          ggtitle("Bias of Treatment Coefficient", subtitle=subtitle) +
          ylab("") +
          xlab(xlab)

  return(plot)
}






















#
