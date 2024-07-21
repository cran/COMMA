#' Predictive Value Weighting Estimation of the Binary Mediator Misclassification Model
#' 
#' Estimate \eqn{\beta}, \eqn{\gamma}, and \eqn{\theta} parameters from 
#' the true mediator, observed mediator, and outcome mechanisms, respectively,
#' in a binary mediator misclassification model using a predictive value weighting
#' approach.
#' 
#' Note that this method can only be used for binary outcome models.
#'
#' @param Mstar A numeric vector of indicator variables (1, 2) for the observed
#'   mediator \code{M*}. There should be no \code{NA} terms. The reference category is 2.
#' @param outcome A vector containing the outcome variables of interest. There
#'   should be no \code{NA} terms.
#' @param outcome_distribution A character string specifying the distribution of 
#'   the outcome variable. Options are \code{"Bernoulli"},
#'   \code{"Poisson"}, or \code{"Normal"}.
#' @param interaction_indicator A logical value indicating if an interaction between
#'   \code{x} and \code{m} should be used to generate the outcome variable, \code{y}.
#' @param x_matrix A numeric matrix of predictors in the true mediator and outcome mechanisms.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param c_matrix A numeric matrix of covariates in the true mediator and outcome mechanisms.
#'   \code{c_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param beta_start A numeric vector or column matrix of starting values for the \eqn{\beta}
#'   parameters in the true mediator mechanism. The number of elements in \code{beta_start}
#'   should be equal to the number of columns of \code{x_matrix} and \code{c_matrix} plus 1.
#'   Starting values should be provided in the following order: intercept, slope 
#'   coefficient for the \code{x_matrix} term, slope coefficient for first column
#'   of the \code{c_matrix}, ..., slope coefficient for the final column of the \code{c_matrix}.
#' @param gamma_start A numeric vector or matrix of starting values for the \eqn{\gamma}
#'   parameters in the observation mechanism. In matrix form, the \code{gamma_start} matrix rows
#'   correspond to parameters for the \code{M* = 1}
#'   observed mediator, with the dimensions of \code{z_matrix} plus 1, and the
#'   gamma parameter matrix columns correspond to the true mediator categories
#'   \eqn{M \in \{1, 2\}}. A numeric vector for \code{gamma_start} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_start <- c(gamma_matrix)}.
#'   Starting values should be provided in the following order within each column:
#'   intercept, slope coefficient for first column of the \code{z_matrix}, ...,
#'   slope coefficient for the final column of the \code{z_matrix}.
#' @param theta_start A numeric vector or column matrix of starting values for the \eqn{\theta}
#'   parameters in the outcome mechanism. The number of elements in \code{theta_start}
#'   should be equal to the number of columns of \code{x_matrix} and \code{c_matrix} plus 2 
#'   (if \code{interaction_indicator} is \code{FALSE}) or 3 (if
#'   \code{interaction_indicator} is \code{TRUE}). Starting values should be
#'   provided in the following order: intercept, slope coefficient for the \code{x_matrix} term,
#'   slope coefficient for the mediator \code{m} term,
#'   slope coefficient for first column of the \code{c_matrix}, ...,
#'   slope coefficient for the final column of the \code{c_matrix},
#'   and, optionally, slope coefficient for \code{xm}).
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param em_method A character string specifying which EM algorithm will be applied.
#'   Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
#'   recommended option is \code{"squarem"}.
#'
#' @return \code{COMMA_PVW} returns a data frame containing four columns. The first
#'   column, \code{Parameter}, represents a unique parameter value for each row.
#'   The next column contains the parameter \code{Estimates}. The third column,
#'   \code{Convergence}, reports whether or not the algorithm converged for a
#'   given parameter estimate. The final column, \code{Method}, reports
#'   that the estimates are obtained from the "PVW" procedure. 
#'   
#' @export
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include COMBO_EM_algorithm.R
#' @include COMBO_EM_function.R
#' @include COMBO_weight.R
#' @include COMMA_data.R
#' 
#' @importFrom stats lm poisson
#'
#' @examples
#' set.seed(20240709)
#' sample_size <- 2000
#' 
#' n_cat <- 2 # Number of categories in the binary mediator
#' 
#' # Data generation settings
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#' c_shape <- 1
#' 
#' # True parameter values (gamma terms set the misclassification rate)
#' true_beta <- matrix(c(1, -2, .5), ncol = 1)
#' true_gamma <- matrix(c(1, 1, -.5, -1.5), nrow = 2, byrow = FALSE)
#' true_theta <- matrix(c(1, 1.5, -2, -.2), ncol = 1)
#' 
#' example_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
#'                            interaction_indicator = FALSE,
#'                            outcome_distribution = "Bernoulli",
#'                            true_beta, true_gamma, true_theta)
#'                            
#' beta_start <- matrix(rep(1, 3), ncol = 1)
#' gamma_start <- matrix(rep(1, 4), nrow = 2, ncol = 2)
#' theta_start <- matrix(rep(1, 4), ncol = 1)
#' 
#' Mstar = example_data[["obs_mediator"]]
#' outcome = example_data[["outcome"]]
#' x_matrix = example_data[["x"]]
#' z_matrix = example_data[["z"]]
#' c_matrix = example_data[["c"]]
#'                            
#' PVW_results <- COMMA_PVW(Mstar, outcome, outcome_distribution = "Bernoulli",
#'                          interaction_indicator = FALSE,
#'                          x_matrix, z_matrix, c_matrix,
#'                          beta_start, gamma_start, theta_start)
#'
#' PVW_results
#' 
COMMA_PVW <- function(Mstar, # Observed mediator vector
                      outcome, # Outcome vector
                      outcome_distribution,
                      interaction_indicator,
                      # Predictor matrices
                      x_matrix, z_matrix, c_matrix,
                      # Start values for parameters
                      beta_start, gamma_start,
                      theta_start,
                      # EM settings
                      tolerance = 1e-7, max_em_iterations = 1500,
                      em_method = "squarem"){
  
  n_cat = 2 # Number of categories in mediator
  sample_size = length(Mstar) # Sample size
  
  # Create design matrices
  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)
  
  x_mat <- as.matrix(x_matrix)
  c_mat <- as.matrix(c_matrix)
  
  # Create matrix of true mediation model predictors
  mediation_model_predictors <- cbind(x_matrix, c_matrix)
  
  # Run the COMBO EM algorithm for the true and observed mediation model
  COMBO_EM_results <- COMBO_EM_algorithm(Mstar,
                                         mediation_model_predictors,
                                         z_matrix,
                                         beta_start, gamma_start,
                                         tolerance, max_em_iterations,
                                         em_method)
  # Save results
  gamma_index_1 = ncol(mediation_model_predictors) + 2
  gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
  
  n_param <- length(c(beta_start, c(gamma_start), theta_start))
  
  predicted_beta <- matrix(COMBO_EM_results$Estimates[1:(ncol(mediation_model_predictors) + 1)],
                           ncol = 1)
  predicted_gamma <- matrix(COMBO_EM_results$Estimates[gamma_index_1:gamma_index_2], 
                            ncol = 2, byrow = FALSE)
  
  # Create a matrix of observed mediator variables using dummy coding
  mstar_matrix <- matrix(c(ifelse(Mstar == 1, 1, 0),
                           ifelse(Mstar == 2, 1, 0)),
                         ncol = 2, byrow = FALSE)
  
  # Create matrix of predictors for the true mediator
  X_design <- cbind(rep(1, sample_size), mediation_model_predictors)
  
  # Generate probabilities for the true mediator value based on EM results
  pi_matrix <- pi_compute(predicted_beta, X_design, sample_size, n_cat)
  
  # Create matrix of predictors for the observed mediator
  Z_design <- matrix(c(rep(1, sample_size), z_matrix),
                     nrow = sample_size, byrow = FALSE)
  
  # Generate probabilities for observed mediator conditional on true mediator
  ## Based on EM results
  pistar_matrix <- pistar_compute(predicted_gamma, Z_design, sample_size, n_cat)
  
  # Estimate sensitivity and specificity
  sensitivity <- pistar_matrix[1:sample_size, 1]
  specificity <- pistar_matrix[(sample_size + 1):(2 * sample_size), 2]
  
  # Organize data for model predicting the observed mediator
  mstar_model_data <- data.frame(x = x_matrix, c = c_matrix, z = z_matrix,
                                 y = outcome,
                                 mstar_01 = ifelse(Mstar == 1, 1, 0))
  
  # Fit model for observed mediator based on x, c, y, z (with interactions)
  mstar_model <- glm(mstar_01 ~ . ^2,
                     data = mstar_model_data, family = "binomial")
  
  # Predict observed mediators
  predictions <- stats::predict(mstar_model, type = "response")
  
  # Ensure no exact 0 or 1 values
  sensitivity[predictions >= sensitivity] <- predictions[predictions >= sensitivity] + 0.001
  specificity[predictions <= (1-specificity)] <- 1 - predictions[predictions <= (1 - specificity)] + 0.001
  
  # Compute NPV and PPV
  term1 <- (sensitivity - 1) * predictions * (1 / (sensitivity * (predictions - 1)))
  term2 <- (specificity - 1) * (predictions - 1) * (1 / (specificity * predictions))
  det <- 1/(term1*term2-1)
  ppv_calc <- det * (term2 - 1)
  npv_calc <- det * (term1 - 1)
  
  ppv <- unname(ppv_calc)
  npv <- unname(npv_calc)
  
  if(interaction_indicator == FALSE & outcome_distribution == "Bernoulli"){
    
    # Duplicate the dataset
    actual_dataset <- data.frame(x = x_matrix, m = 0, c = c_matrix, 
                                 y_01 = ifelse(outcome == 1, 1, 0),
                                 mstar_01 = mstar_model_data$mstar_01)
    
    duplicate_dataset <- data.frame(x = x_matrix, m = 1, c = c_matrix,
                                    y_01 = ifelse(outcome == 1, 1, 0),
                                    mstar_01 = mstar_model_data$mstar_01)
    
    doubled_data <- rbind(actual_dataset, duplicate_dataset)
    
    # Apply NPV and PPV weights
    doubled_data$w <- 0
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 1] <- ppv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 1) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 1] <- 1 - ppv[doubled_data$m == 0 & doubled_data$mstar_01 == 1]
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 0] <- 1 - npv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 0) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 0] <- npv[doubled_data$m == 0 & doubled_data$mstar_01 == 0]
    
    # Remove mstar term from dataset before modeling.
    doubled_data_2 <- doubled_data[,-(ncol(doubled_data) - 1)]
    
    # Fit weighted logistic regression to estimate theta
    weighted_outcome_model <- glm(y_01 ~ . -y_01 -w, weights = w,
                                  data = doubled_data_2,
                                  family = "binomial"(link = "logit"))
    summary(weighted_outcome_model)
    
    #  Save results
    beta_param_names <- paste0(rep("beta_", ncol(X_design)), 1:ncol(X_design))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)))
    
    
    PVW_results <- data.frame(Parameter = c(beta_param_names,
                                            gamma_param_names,
                                            theta_param_names),
                              Estimates = c(c(predicted_beta),
                                            c(predicted_gamma),
                                            c(unname(coefficients(weighted_outcome_model)))),
                              Convergence = rep(COMBO_EM_results$Convergence[1],
                                                n_param),
                              Method = "PVW")
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Bernoulli"){
    
    # Duplicate the dataset
    interaction_m0 = x_matrix * 0
    actual_dataset <- data.frame(x = x_matrix, m = 0, c = c_matrix,
                                 xm = interaction_m0,
                                 y_01 = ifelse(outcome == 1, 1, 0),
                                 mstar_01 = mstar_model_data$mstar_01)
    
    interaction_m1 = x_matrix * 1
    duplicate_dataset <- data.frame(x = x_matrix, m = 1, c = c_matrix,
                                    xm = interaction_m1,
                                    y_01 = ifelse(outcome == 1, 1, 0),
                                    mstar_01 = mstar_model_data$mstar_01)
    
    doubled_data <- rbind(actual_dataset, duplicate_dataset)
    
    # Apply NPV and PPV weights
    doubled_data$w <- 0
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 1] <- ppv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 1) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 1] <- 1 - ppv[doubled_data$m == 0 & doubled_data$mstar_01 == 1]
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 0] <- 1 - npv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 0) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 0] <- npv[doubled_data$m == 0 & doubled_data$mstar_01 == 0]
    
    # Remove mstar term from dataset before modeling.
    doubled_data_2 <- doubled_data[,-(ncol(doubled_data) - 1)]
    
    w <- doubled_data_2$w
    
    # Fit weighted logistic regression to estimate theta
    weighted_outcome_model <- glm(y_01 ~ . -y_01 -w, weights = w,
                                  data = doubled_data_2,
                                  family = "binomial"(link = "logit"))
    summary(weighted_outcome_model)
    
    #  Save results
    beta_param_names <- paste0(rep("beta_", ncol(X_design)), 1:ncol(X_design))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)),
                           "theta_xm")
    
    
    PVW_results <- data.frame(Parameter = c(beta_param_names,
                                            gamma_param_names,
                                            theta_param_names),
                              Estimates = c(c(predicted_beta),
                                            c(predicted_gamma),
                                            c(unname(coefficients(weighted_outcome_model)))),
                              Convergence = rep(COMBO_EM_results$Convergence[1],
                                                n_param),
                              Method = "PVW")
    
  } else if(interaction_indicator == FALSE & outcome_distribution == "Poisson"){
    
    # Duplicate the dataset
    actual_dataset <- data.frame(x = x_matrix, m = 0, c = c_matrix, 
                                 y = outcome,
                                 mstar_01 = mstar_model_data$mstar_01)
    
    duplicate_dataset <- data.frame(x = x_matrix, m = 1, c = c_matrix,
                                    y = outcome,
                                    mstar_01 = mstar_model_data$mstar_01)
    
    doubled_data <- rbind(actual_dataset, duplicate_dataset)
    
    # Apply NPV and PPV weights
    doubled_data$w <- 0
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 1] <- ppv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 1) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 1] <- 1 - ppv[doubled_data$m == 0 & doubled_data$mstar_01 == 1]
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 0] <- 1 - npv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 0) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 0] <- npv[doubled_data$m == 0 & doubled_data$mstar_01 == 0]
    
    # Remove mstar term from dataset before modeling.
    doubled_data_2 <- doubled_data[,-(ncol(doubled_data) - 1)]
    
    # Remove negative weights (why are these happening?)
    doubled_data_2$w_no_negative <- ifelse(doubled_data_2$w < 0, 0, doubled_data_2$w)
    
    w_no_negative <- doubled_data_2$w_no_negative
    
    # Fit weighted logistic regression to estimate theta
    weighted_outcome_model <- glm(y ~ . -y -w -w_no_negative, weights = w_no_negative,
                                  data = doubled_data_2,
                                  family = "poisson"(link = "log"))
    summary(weighted_outcome_model)
    
    #  Save results
    beta_param_names <- paste0(rep("beta_", ncol(X_design)), 1:ncol(X_design))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)))
    
    
    PVW_results <- data.frame(Parameter = c(beta_param_names,
                                            gamma_param_names,
                                            theta_param_names),
                              Estimates = c(c(predicted_beta),
                                            c(predicted_gamma),
                                            c(unname(coefficients(weighted_outcome_model)))),
                              Convergence = rep(COMBO_EM_results$Convergence[1],
                                                n_param),
                              Method = "PVW")
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Poisson"){
    
    # Duplicate the dataset
    interaction_m0 = x_matrix * 0
    actual_dataset <- data.frame(x = x_matrix, m = 0, c = c_matrix,
                                 xm = interaction_m0,
                                 y = outcome,
                                 mstar_01 = mstar_model_data$mstar_01)
    
    interaction_m1 = x_matrix * 1
    duplicate_dataset <- data.frame(x = x_matrix, m = 1, c = c_matrix,
                                    xm = interaction_m1,
                                    y = outcome,
                                    mstar_01 = mstar_model_data$mstar_01)
    
    doubled_data <- rbind(actual_dataset, duplicate_dataset)
    
    # Apply NPV and PPV weights
    doubled_data$w <- 0
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 1] <- ppv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 1) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 1] <- 1 - ppv[doubled_data$m == 0 & doubled_data$mstar_01 == 1]
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 0] <- 1 - npv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 0) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 0] <- npv[doubled_data$m == 0 & doubled_data$mstar_01 == 0]
    
    # Remove mstar term from dataset before modeling.
    doubled_data_2 <- doubled_data[,-(ncol(doubled_data) - 1)]
    
    w <- doubled_data_2$w
    
    # Remove negative weights (why are these happening?)
    doubled_data_2$w_no_negative <- ifelse(doubled_data_2$w < 0, 0, doubled_data_2$w)
    
    w_no_negative <- doubled_data_2$w_no_negative
    
    # Fit weighted logistic regression to estimate theta
    weighted_outcome_model <- glm(y ~ . -y -w -w_no_negative, weights = w_no_negative,
                                  data = doubled_data_2,
                                  family = "poisson"(link = "log"))
    summary(weighted_outcome_model)
    
    #  Save results
    beta_param_names <- paste0(rep("beta_", ncol(X_design)), 1:ncol(X_design))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)),
                           "theta_xm")
    
    
    PVW_results <- data.frame(Parameter = c(beta_param_names,
                                            gamma_param_names,
                                            theta_param_names),
                              Estimates = c(c(predicted_beta),
                                            c(predicted_gamma),
                                            c(unname(coefficients(weighted_outcome_model)))),
                              Convergence = rep(COMBO_EM_results$Convergence[1],
                                                n_param),
                              Method = "PVW")
    
  } else if(interaction_indicator == FALSE & outcome_distribution == "Normal"){
    
    # Duplicate the dataset
    actual_dataset <- data.frame(x = x_matrix, m = 0, c = c_matrix, 
                                 y = outcome,
                                 mstar_01 = mstar_model_data$mstar_01)
    
    duplicate_dataset <- data.frame(x = x_matrix, m = 1, c = c_matrix,
                                    y = outcome,
                                    mstar_01 = mstar_model_data$mstar_01)
    
    doubled_data <- rbind(actual_dataset, duplicate_dataset)
    
    # Apply NPV and PPV weights
    doubled_data$w <- 0
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 1] <- ppv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 1) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 1] <- 1 - ppv[doubled_data$m == 0 & doubled_data$mstar_01 == 1]
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 0] <- 1 - npv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 0) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 0] <- npv[doubled_data$m == 0 & doubled_data$mstar_01 == 0]
    
    # Remove mstar term from dataset before modeling.
    doubled_data_2 <- doubled_data[,-(ncol(doubled_data) - 1)]
    
    # Remove negative weights (why are these happening?)
    doubled_data_2$w_no_negative <- ifelse(doubled_data_2$w < 0, 0, doubled_data_2$w)
    
    w_no_negative <- doubled_data_2$w_no_negative
    
    # Fit weighted regression to estimate theta
    weighted_outcome_model <- lm(y ~ . -y -w -w_no_negative, weights = w_no_negative,
                                 data = doubled_data_2)
    summary(weighted_outcome_model)
    
    
    #  Save results
    beta_param_names <- paste0(rep("beta_", ncol(X_design)), 1:ncol(X_design))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)))
    
    
    PVW_results <- data.frame(Parameter = c(beta_param_names,
                                            gamma_param_names,
                                            theta_param_names),
                              Estimates = c(c(predicted_beta),
                                            c(predicted_gamma),
                                            c(unname(coefficients(weighted_outcome_model)))),
                              Convergence = rep(COMBO_EM_results$Convergence[1],
                                                n_param),
                              Method = "PVW")
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Normal"){
    
    # Duplicate the dataset
    interaction_m0 = x_matrix * 0
    actual_dataset <- data.frame(x = x_matrix, m = 0, c = c_matrix,
                                 xm = interaction_m0,
                                 y = outcome,
                                 mstar_01 = mstar_model_data$mstar_01)
    
    interaction_m1 = x_matrix * 1
    duplicate_dataset <- data.frame(x = x_matrix, m = 1, c = c_matrix,
                                    xm = interaction_m1,
                                    y = outcome,
                                    mstar_01 = mstar_model_data$mstar_01)
    
    doubled_data <- rbind(actual_dataset, duplicate_dataset)
    
    # Apply NPV and PPV weights
    doubled_data$w <- 0
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 1] <- ppv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 1) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 1] <- 1 - ppv[doubled_data$m == 0 & doubled_data$mstar_01 == 1]
    doubled_data$w[doubled_data$m == 1 & doubled_data$mstar_01 == 0] <- 1 - npv[which(doubled_data$m == 1 & doubled_data$mstar_01 == 0) - sample_size]
    doubled_data$w[doubled_data$m == 0 & doubled_data$mstar_01 == 0] <- npv[doubled_data$m == 0 & doubled_data$mstar_01 == 0]
    
    # Remove mstar term from dataset before modeling.
    doubled_data_2 <- doubled_data[,-(ncol(doubled_data) - 1)]
    
    w <- doubled_data_2$w
    
    # Remove negative weights (why are these happening?)
    doubled_data_2$w_no_negative <- ifelse(doubled_data_2$w < 0, 0, doubled_data_2$w)
    
    w_no_negative <- doubled_data_2$w_no_negative
    
    # Fit weighted regression to estimate theta
    weighted_outcome_model <- lm(y ~ . -y -w -w_no_negative, weights = w_no_negative,
                                  data = doubled_data_2)
    summary(weighted_outcome_model)
    
    #  Save results
    beta_param_names <- paste0(rep("beta_", ncol(X_design)), 1:ncol(X_design))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)),
                           "theta_xm")
    
    
    PVW_results <- data.frame(Parameter = c(beta_param_names,
                                            gamma_param_names,
                                            theta_param_names),
                              Estimates = c(c(predicted_beta),
                                            c(predicted_gamma),
                                            c(unname(coefficients(weighted_outcome_model)))),
                              Convergence = rep(COMBO_EM_results$Convergence[1],
                                                n_param),
                              Method = "PVW")
    
  } else {
    
    "Undefined"
  }
  
  return(PVW_results)
}
