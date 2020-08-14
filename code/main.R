### CompStat project --- Max Schäfer

rm(list=ls())


#install.packages("glmnet")
#install.packages("randomcoloR")

library(glmnet)
library(mvtnorm)

library(stringr)

library(ggplot2)
library(randomcoloR)

# Own modules
#source("estimation.R")
#source("data_generating_process.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/estimation.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/data_generating_process.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/variable_selection.R")
source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/compute_metrics.R")




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


source("C:/Users/Max/Desktop/MASTER/computational_statistics/compstat-project/code/data_generating_process.R")

### Main ---
set.seed(1234)

R <- 50

selection_method <- "double"

n <- 200
n_F_attr <- 30
n_G_attr <- 50
n_H_attr <- 30

beta_F_size <- .25
beta_H_size <- 1
beta_GD_size <- 1
beta_GY_size <- .5
treatment_effect <- .5

unconfoundedness_rate <- .9


# Special case
n <- 200
n_G_attr <- 100
beta_DY <- .5



# Housekeeping of feature names and identifiers for identification.
# Create vector with covariate names.
colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))

# Required for first selection step
colnames_DG <- c("D", colname_G)
# Required for second selection step
colnames_GH <- c(colname_G, colname_H)
colnumbers_G <- seq(from=1, to=n_G_attr, by=1)
colnumbers_GH <- seq(from=1, to=(n_G_attr+n_H_attr), by=1)



# Main : a single simulation run
sim_results_vec <- replicate(R, compute_metrics(n=n,
                                     n_F_attr=n_F_attr,
                                     n_G_attr=n_G_attr,
                                     n_H_attr=n_H_attr,
                                     beta_GD_size=beta_GD_size,
                                     beta_GY_size=beta_GY_size,
                                     beta_H_size=beta_H_size,
                                     beta_F_size=beta_F_size,
                                     treatment_effect=treatment_effect,
                                     unconfoundedness_rate=unconfoundedness_rate,
                                     selection_method=selection_method))

confounder_bias_vec <- unlist(sim_results_vec[1, 1:R])
test_MSE_vec <- unlist(sim_results_vec[2, 1:R])

tp_selection_rate_G <- unlist(sim_results_vec[3, 1:R])
tn_selection_rate_G <- unlist(sim_results_vec[4, 1:R])
fp_selection_rate_G <- unlist(sim_results_vec[5, 1:R])
fn_selection_rate_G <- unlist(sim_results_vec[6, 1:R])

select_treatment_identifier <- unlist(sim_results_vec[7, 1:R]) #TODO

# Evaluate
cat("Mean test MSE is", mean(test_MSE_vec))
cat("Median confounder bias is", median(confounder_bias_vec))
cat("Mean of false positive 'predictions' is", mean(fp_selection_rate_G), "i.e. the more false positives the larger the confounding problem")
hist(confounder_bias_vec)
cat("RMSE is", sqrt( mean( confounder_bias_vec**2 )))

df <- data.frame(confounder_bias_vec)
plot <- ggplot(data=df, aes_string(x="confounder_bias_vec")) +
        geom_histogram(aes(y=..density..), binwidth=.02, colour="blue", fill="white") +
        geom_density(alpha=.2, fill="#FF6666")
plot


# Simulations

# Vary over
iter_over <- "unconfoundedness_rate"

unconfoundedness_rate_vec <- seq(from=.8, to=1, by=.025)

treatment_effect_vec <- seq(from=0, to=3, by=.2)
beta_GD_size_vec <- seq(from=0, to=4, by=.2)
beta_GY_size_vec <- seq(from=0, to=2, by=.1)
beta_H_size_vec <- seq(from=0, to=2, by=.1)
beta_F_size_vec <- seq(from=0, to=.5, by=.05)

n_G_attr_vec <- seq(from=10, to=120, by=10)
n_F_attr_vec <- seq(from=10, to=50, by=5)

if (iter_over=="unconfoundedness_rate") {
  sim_parameter_vec <- unconfoundedness_rate_vec
}
if (iter_over=="treatment_effect") {
  sim_parameter_vec <- treatment_effect_vec
}
if (iter_over=="n_G_attr") {
  sim_parameter_vec <- n_G_attr_vec
}
if (iter_over=="n_F_attr") {
  sim_parameter_vec <- n_F_attr_vec
}
if (iter_over=="beta_F_size") {
  sim_parameter_vec <- beta_F_size_vec
}
if (iter_over=="beta_H_size") {
    sim_parameter_vec <- beta_H_size_vec
}
if (iter_over=="beta_GY_size") {
  sim_parameter_vec <- beta_GY_size_vec
}



# Collect simulation results
simulation_results <- simulation_wrapper(R=R, iter_over=iter_over,
                              sim_parameter_vec=sim_parameter_vec,
                              selection_method=selection_method,
                              n=n, n_F_attr=n_F_attr,
                              n_G_attr=n_G_attr, n_H_attr=n_H_attr,
                              treatment_effect=treatment_effect,
                              beta_GD_size=beta_GD_size,
                              beta_GY_size=beta_GY_size,
                              beta_F_size=beta_F_size,
                              unconfoundedness_rate=unconfoundedness_rate)



mean_squared_bias_vec <- simulation_results$mean_squared_bias_vec
mean_abs_dev_bias_vec <- simulation_results$mean_abs_dev_bias_vec
median_bias_vec <- simulation_results$median_bias_vec
mean_test_MSE_vec <- simulation_results$mean_test_MSE_vec
tp_selection_rate_G_vec <- simulation_results$tp_selection_rate_G_vec
tn_selection_rate_G_vec <- simulation_results$tn_selection_rate_G_vec
fp_selection_rate_G_vec <- simulation_results$fp_selection_rate_G_vec
fn_selection_rate_G_vec <- simulation_results$fn_selection_rate_G_vec
select_treatment_identifier_vec <- simulation_results$select_treatment_identifier_vec


### Visualization
# Metrics --- confounder bias
x_axis <- sim_parameter_vec
df <- data.frame(cbind(x_axis, mean_squared_bias_vec, mean_abs_dev_bias_vec, median_bias_vec))

# Plot
plot <- ggplot(data=df, aes(x=x_axis))
plot <- plot +
        geom_line(aes_string(y="mean_squared_bias_vec"), size=1, col="blue", alpha=1) +
        geom_line(aes_string(y="mean_abs_dev_bias_vec"), size=1, col="red", alpha=1) +
        geom_line(aes_string(y="median_bias_vec"), size=1, col="grey", alpha=1) +
        geom_line(aes(y=0), lty="dashed") +
        ylab("Bias of Treatment Coefficient") +
        xlab(iter_over) #n_G_attr

plot


# Metrics --- prediction
x_axis <- sim_parameter_vec
df <- data.frame(cbind(x_axis, mean_test_MSE_vec))

# Plot
plot <- ggplot(data=df, aes(x=x_axis))
plot <- plot +
        geom_line(aes_string(y="mean_test_MSE_vec"), size=1, col="blue", alpha=1) +
        geom_line(aes(y=0), lty="dashed") +
        ylab("Mean test MSE (prediction)") +
        xlab(iter_over)

plot


# Metrics --- variable selection
x_axis <- sim_parameter_vec
df <- data.frame(cbind(x_axis,
                       tp_selection_rate_G_vec, tn_selection_rate_G_vec,
                       fp_selection_rate_G_vec, fn_selection_rate_G_vec))

# Plot
plot <- ggplot(data=df, aes(x=x_axis))
plot <- plot +
        geom_line(aes_string(y="tp_selection_rate_G_vec"), size=1, col="blue", alpha=.9) +
        geom_line(aes_string(y="tn_selection_rate_G_vec"), size=1, col="blue", alpha=.5) +
        geom_line(aes_string(y="fp_selection_rate_G_vec"), size=1, col="red", alpha=.9) +
        geom_line(aes_string(y="fn_selection_rate_G_vec"), size=1, col="red", alpha=.5) +
        geom_line(aes(y=0), lty="dashed") +
        ylab("TP, TN, FP, FN rates for selecting confounder") +
        xlab(iter_over)

plot




library(tidyr)

df_long <- tidyr::gather(df[ , 2:dim(df)[2]], key="metric")
df_long["x_axis"] <- df[,1]
plot <- ggplot(data=df_long, aes(x=x_axis, y=value, colour=metric)) +
        geom_line() +
        geom_line(aes(y=0), lty="dashed") +
        ylab("TP, TN, FP, FN rates for selecting confounder") +
        xlab("Unconfoundedness Rate")

plot



# 2. Use Lasso to select covariates given their association with treatment D
#TODO cannot use logistic model then?? or use extension of lasso with logit...
# type.measure class gives misclassification error for cv
#lasso_two_cv <- cv.glmnet(G_train, D_train, alpha=1, family="binomial", type.measure="class", intercept=FALSE)














#TODO Quick inspection
#plot(lasso_one_cv)
#plot(lasso_one, label=TRUE)
#print(lasso_one)
#coef(lasso_one, s=.1)

# >>>
# Quick look
n_F_attr <- 12
n_G_attr <- 12
n_H_attr <- 12
treatment_effect <- 1
unconfoundedness_rate <- .9


colname_F <- str_c(rep("F_", n_F_attr), seq(from=1, to=n_F_attr, by=1))
colname_G <- str_c(rep("G_", n_G_attr), seq(from=1, to=n_G_attr, by=1))
colname_H <- str_c(rep("H_", n_H_attr), seq(from=1, to=n_H_attr, by=1))


data_set <- generate_data(n=n,
                          n_F_attr=n_F_attr,
                          n_G_attr=n_G_attr,
                          n_H_attr=n_H_attr,
                          treatment_effect=treatment_effect,
                          beta_GY_size = 2,
                          unconfoundedness_rate=unconfoundedness_rate)

data_train <- data_set$data
y_train <- data_train$y
D_train <- data_train$D
F_train <- data.matrix(data_train[colname_F])
G_train <- data.matrix(data_train[colname_G])
H_train <- data.matrix(data_train[colname_H])


# Assemble different regressor matrices
#X_train <- data.matrix(cbind(D_train, F_train, G_train, H_train))
X_train <- data.matrix(cbind(D_train, G_train, H_train))
Z_train <- data.matrix(cbind(G_train, H_train))
#Z_train <- data.matrix(cbind(F_train, G_train, H_train))

# Regression results
#lm(y_train ~ D_train + F_train + G_train + H_train)
lm(y_train ~ D_train + G_train + H_train)
lm(y_train ~ G_train + H_train)
lm(y_train ~ D_train + H_train)

lm(y_train ~ D_train + G_train[, 1:12] + H_train)
lm(y_train ~ D_train)

#lm(D_train ~ F_train + G_train + H_train)
lm(D_train ~ G_train + H_train)
lm(D_train ~ G_train)

#beta_hat <- least_squares_estimator(y=y_train, X=X_train)
# <<<




















### Code snippets

# Ensure correlations are positive
#A <- matrix(rnorm(n=H_no_attr**2, mean=0, sd=2), ncol=H_no_attr)
#sigma <- t(A) %*% A
#diag(sigma) <- diag(sigma) * 10
#sigma <- abs(sigma)

#x <- rbinom(n=10, size=100, prob=.15)
#print(x)


# Initialize plot.
#plot <- ggplot(data=df, aes(x=lambda_grid))
# Loop through each beta and plot its value with increasing penalty.
#for (name in cols) {
#  plot <- plot +
#          geom_line(aes_string(y=name), size=1, col=randomColor(count=1, luminosity="light"))
#}
#plot <- plot +
#        geom_line(aes(y=0), lty="dashed") +
#        ylab("Lasso coefficients") +
#        xlab("Lambda")
#plot




#
