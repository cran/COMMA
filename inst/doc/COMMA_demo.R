## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE, error = FALSE, message = FALSE,
                      fig.align = "center")

library(ggplot2)
library(kableExtra)

#devtools::install_github("kimberlywebb/COMMA")

## -----------------------------------------------------------------------------
library(COMMA)
library(dplyr)

set.seed(20240422)

sample_size <- 10000

n_cat <- 2 # Number of categories in the binary mediator

# Data generation settings
x_mu <- 0
x_sigma <- 1
z_shape <- 1
z_scale <- 1
c_shape <- 1

# True parameter values (gamma terms set the misclassification rate)
true_beta <- matrix(c(1, -2, .5), ncol = 1)
true_gamma <- matrix(c(1.8, 1, -1.5, -1), nrow = 2, byrow = FALSE) 
true_theta <- matrix(c(1, 1.5, -2.5, -.2), ncol = 1)

# Generate data.
my_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
                      interaction_indicator = FALSE,
                      outcome_distribution = "Normal",
                      true_beta, true_gamma, true_theta)

# Save list elements as vectors.
Mstar = my_data[["obs_mediator"]]
Mstar_01 <- ifelse(Mstar == 1, 1, 0)
outcome = my_data[["outcome"]]
x_matrix = my_data[["x"]]
z_matrix = my_data[["z"]]
c_matrix = my_data[["c"]]

## -----------------------------------------------------------------------------
# Supply starting values for all parameters.
beta_start <- coef(glm(Mstar_01 ~ x_matrix + c_matrix,
                       family = "binomial"(link = "logit")))
gamma_start <- matrix(rep(1,4), ncol = 2, nrow = 2, byrow = FALSE)
theta_start <- coef(lm(outcome ~ x_matrix + Mstar_01 + c_matrix))

# Estimate parameters using the EM-Algorithm.
EM_results <- COMMA_EM(Mstar, outcome, outcome_distribution = "Normal",
                       interaction_indicator = FALSE,
                       x_matrix, z_matrix, c_matrix,
                       beta_start, gamma_start, theta_start, sigma_start = 1)

EM_results$True_Value <- c(true_beta, c(true_gamma), true_theta, 1)
EM_results$Estimates <- round(EM_results$Estimates, 3)
EM_results

## -----------------------------------------------------------------------------
# Estimate parameters using the OLS correction.
OLS_results <- COMMA_OLS(Mstar, outcome,
                         x_matrix, z_matrix, c_matrix,
                         beta_start, gamma_start, theta_start)

OLS_results$True_Value <- c(true_beta, c(true_gamma), true_theta[c(1,3,2,4)])
OLS_results$Estimates <- round(OLS_results$Estimates, 3)
OLS_results

## -----------------------------------------------------------------------------
NormalY_results <- data.frame(Parameter = rep(c("beta_x", "theta_x", "theta_m"),
                                              3),
                              Method = c(rep("EM", 3), rep("OLS", 3),
                                         rep("Naive", 3)),
                              Estimate = c(EM_results$Estimates[c(2, 9, 10)],
                                           OLS_results$Estimates[c(2, 9, 10)],
                                           beta_start[2], theta_start[c(2,3)]),
                              True_Value = rep(c(true_beta[2],
                                                 true_theta[2], true_theta[3]),
                                               3))

ggplot(data = NormalY_results) +
  geom_point(aes(x = Method, y = Estimate, color = Method),
             size = 3) +
  geom_hline(aes(yintercept = True_Value),
             linetype = "dashed") +
  facet_wrap(~Parameter) +
  theme_minimal() +
  scale_color_manual(values = c("#DA86A5", "#785E95", "#409DBE")) +
  ggtitle("Parameter estimates from EM, OLS, and Naive approaches",
          subtitle = "Normal outcome model with a misclassified mediator")

## -----------------------------------------------------------------------------
library(COMMA)
library(dplyr)

set.seed(20240423)

sample_size <- 10000

n_cat <- 2 # Number of categories in the binary mediator

# Data generation settings
x_mu <- 0
x_sigma <- 1
z_shape <- 1
z_scale <- 1
c_shape <- 1

# True parameter values (gamma terms set the misclassification rate)
true_beta <- matrix(c(1, -2, .5), ncol = 1)
true_gamma <- matrix(c(1.8, 1, -1.5, -1), nrow = 2, byrow = FALSE) 
true_theta <- matrix(c(1, 1.5, -2.5, -.2, .5), ncol = 1)

# Generate data.
my_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
                      interaction_indicator = TRUE,
                      outcome_distribution = "Bernoulli",
                      true_beta, true_gamma, true_theta)

# Save list elements as vectors.
Mstar = my_data[["obs_mediator"]]
Mstar_01 <- ifelse(Mstar == 1, 1, 0)
outcome = my_data[["outcome"]]
x_matrix = my_data[["x"]]
z_matrix = my_data[["z"]]
c_matrix = my_data[["c"]]

## -----------------------------------------------------------------------------
# Supply starting values for all parameters.
beta_start <- coef(glm(Mstar_01 ~ x_matrix + c_matrix,
                       family = "binomial"(link = "logit")))
gamma_start <- matrix(rep(1,4), ncol = 2, nrow = 2, byrow = FALSE)

xm_interaction <- x_matrix * c_matrix
theta_start <- coef(glm(outcome ~ x_matrix + Mstar_01 + c_matrix +
                          xm_interaction,
                        family = "binomial"(link = "logit")))

# Estimate parameters using the EM-Algorithm.
EM_results <- COMMA_EM(Mstar, outcome, outcome_distribution = "Bernoulli",
                       interaction_indicator = TRUE,
                       x_matrix, z_matrix, c_matrix,
                       beta_start, gamma_start, theta_start)

EM_results$True_Value <- c(true_beta, c(true_gamma), true_theta)
EM_results$Estimates <- round(EM_results$Estimates, 3)
EM_results

## -----------------------------------------------------------------------------
PVW_results <- COMMA_PVW(Mstar, outcome, outcome_distribution = "Bernoulli",
                         interaction_indicator = TRUE,
                         x_matrix, z_matrix, c_matrix,
                         beta_start, gamma_start, theta_start)

PVW_results$True_Value <- c(true_beta, c(true_gamma), true_theta)
PVW_results$Estimates <- round(PVW_results$Estimates, 3)
PVW_results

## -----------------------------------------------------------------------------
BernoulliY_results <- data.frame(Parameter = rep(c("beta_x", "theta_x",
                                                   "theta_m", "theta_xm"),
                                                 3),
                                 Method = c(rep("EM", 4), rep("PVW", 4),
                                            rep("Naive", 4)),
                                 Estimate = c(EM_results$Estimates[c(2, 9, 10, 12)],
                                              PVW_results$Estimates[c(2, 9, 10, 12)],
                                              beta_start[2],
                                              theta_start[c(2,3,5)]),
                                 True_Value = rep(c(true_beta[2],
                                                    true_theta[2],
                                                    true_theta[3],
                                                    true_theta[5]),
                                                  3))

ggplot(data = BernoulliY_results) +
  geom_point(aes(x = Method, y = Estimate, color = Method),
             size = 3) +
  geom_hline(aes(yintercept = True_Value),
             linetype = "dashed") +
  facet_wrap(~Parameter) +
  theme_minimal() +
  scale_color_manual(values = c("#DA86A5", "#785E95", "#ECA698")) +
  ggtitle("Parameter estimates from EM, PVW, and Naive approaches",
          subtitle = "Bernoulli outcome model with a misclassified mediator")

