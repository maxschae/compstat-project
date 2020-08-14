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
