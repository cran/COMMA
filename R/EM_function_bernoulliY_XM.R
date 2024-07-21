#' EM Algorithm Function for Estimation of the Misclassification Model
#' 
#' Function is for cases with \eqn{Y \sim Bernoulli} and with an interaction term
#' in the outcome mechanism.
#'
#' @param param_current A numeric vector of regression parameters, in the order
#'   \eqn{\beta, \gamma, \theta}. The \eqn{\gamma} vector is obtained from the matrix form.
#'   In matrix form, the gamma parameter matrix rows
#'   correspond to parameters for the \code{M* = 1}
#'   observed mediator, with the dimensions of \code{Z}.
#'   In matrix form, the gamma parameter matrix columns correspond to the true mediator categories
#'   \eqn{j = 1, \dots,} \code{n_cat}. The numeric vector \code{gamma_v} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_v <- c(gamma_matrix)}. 
#' @param obs_mediator A numeric vector of indicator variables (1, 2) for the observed
#'   mediator \code{M*}. There should be no \code{NA} terms. The reference category is 2.
#' @param obs_outcome A vector containing the outcome variables of interest. There
#'   should be no \code{NA} terms.
#' @param X A numeric design matrix for the true mediator mechanism.
#' @param Z A numeric design matrix for the observation mechanism.
#' @param c_matrix A numeric matrix of covariates in the true mediator and outcome mechanisms.
#'   \code{c_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param sample_size An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{X} or \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{M},
#'   and the observed outcome, \code{M*} can take.
#'
#' @return \code{EM_function_bernoulliY} returns a numeric vector of updated parameter
#'   estimates from one iteration of the EM-algorithm.
#'   
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include w_m_binaryY.R
#' @include w_m_normalY.R
#' 
#' @importFrom stats rnorm rgamma rmultinom coefficients binomial
#'
EM_function_bernoulliY_XM <- function(param_current,
                                      obs_mediator, obs_outcome,
                                      X, Z, c_matrix,
                                      sample_size, n_cat){
  
  # Create design matrix for true mediator model
  design_matrix = cbind(X, c_matrix)
  
  # Set up parameter indices
  gamma_index_1 = ncol(design_matrix) + 1
  gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
  
  n_param <- length(param_current)

  # Separate current parameters into beta, gamma, theta, sigma vectors
  beta_current = matrix(param_current[1:ncol(design_matrix)], ncol = 1)
  gamma_current = matrix(c(param_current[gamma_index_1:gamma_index_2]),
                         ncol = n_cat, byrow = FALSE)
  theta_current = matrix(c(param_current[(gamma_index_2 + 1):n_param]),
                         ncol = 1)
  
  # Compute probability of each latent mediator value
  probabilities = pi_compute(beta_current, design_matrix, sample_size, n_cat)
  
  # Compute probability of observed mediator, given latent mediator
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)
  
  # Compute likelihood value of Y based on x, m, c, theta, and sigma
  interaction_term_m0 <- X[,-1] * 0
  outcome_design_matrix_m0 <- cbind(cbind(X, cbind(rep(0, sample_size), c_matrix)),
                                    interaction_term_m0)
  model_y_m0 <- outcome_design_matrix_m0 %*% theta_current
  exp_model_y_m0 = exp(model_y_m0)
  p_result_m0 = exp_model_y_m0 / (1 + exp_model_y_m0)
  p_yi_m0 = matrix(c(p_result_m0, 1 - p_result_m0),
                   ncol = 2, nrow = sample_size,
                   byrow = FALSE)
  
  interaction_term_m1 <- X[,-1] * 1
  outcome_design_matrix_m1 <- cbind(cbind(X, cbind(rep(1, sample_size), c_matrix)),
                                    interaction_term_m1)
  model_y_m1 <- outcome_design_matrix_m1 %*% theta_current
  exp_model_y_m1 = exp(model_y_m1)
  p_result_m1 = exp_model_y_m1 / (1 + exp_model_y_m1)
  p_yi_m1 = matrix(c(p_result_m1, 1 - p_result_m1),
                   ncol = 2, nrow = sample_size,
                   byrow = FALSE)

  # Create a matrix of observed mediator variables using dummy coding
  mstar_matrix = matrix(c(ifelse(obs_mediator == 1, 1, 0), 
                          ifelse(obs_mediator == 2, 1, 0)),
                        nrow = sample_size, byrow = FALSE)
  
  # Create a matrix of outcomes
  outcome_matrix = matrix(c(obs_outcome,
                            1 - obs_outcome),
                          nrow = sample_size, byrow = FALSE)
  
  # Compute E-Step weights
  weights = w_m_binaryY(mstar_matrix, outcome_matrix,
                        pistar_matrix = conditional_probabilities,
                        pi_matrix = probabilities,
                        p_yi_m0, p_yi_m1,
                        sample_size, n_cat)

  # Estimate gamma parameters using weighted logistic regression
  ## Weights from E-Step (split by value of latent mediator, m)
  ## Outcome is the observed mediator
  Mstar01 = mstar_matrix[,1]
  fit.gamma1 <- tryCatch({suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,1],
                           family = "binomial"(link = "logit")) )},
                         error = function(e){
                           suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                                                        #weights = weights[,1],
                                                        family = "binomial"(link = "logit")) )
                         })
  gamma1_new <- tryCatch({unname(coefficients(fit.gamma1))},
                         error = function(e){
                           gamma_current[,1]
                         })

  fit.gamma2 <- tryCatch({suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                                                       weights = weights[,2],
                                                       family = "binomial"(link = "logit")) )},
                         error = function(e){
                           suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                                                        #weights = weights[,2],
                                                        family = "binomial"(link = "logit")) )
                         })
  gamma2_new <- tryCatch({unname(coefficients(fit.gamma2))},
                         error = function(e){
                           gamma_current[,2]
                         })

  # Estimate beta parameters using logistic regression
  ## Outcome is the E-Step weight
  fit.beta <- tryCatch({suppressWarnings( stats::glm(weights[,1] ~ . + 0, as.data.frame(design_matrix),
                         family = stats::binomial()) ) },
                       error = function(e){
                         suppressWarnings( stats::glm(Mstar_01 ~ . + 0,
                                                      as.data.frame(design_matrix),
                                                      family = stats::binomial()) )
                       })
  beta_new <- tryCatch({unname(coefficients(fit.beta))},
                       error = function(e){beta_current})

  # Estimate theta parameters using a weighted logistic regression
  ## Duplicate the data, half has m = 0 and half has m = 1
  ## Weights from E-Step (split by value of latent mediator, m)
  ## Outcome is y
  x_vector = X[,-1]
  data1 = data.frame(x = x_vector, m = 0, c = c_matrix,
                     xm = x_vector * 0,
                     w = weights[,2], y = obs_outcome)
  data2 = data.frame(x = x_vector, m = 1, c = c_matrix,
                     xm = x_vector * 1,
                     w = weights[,1], y = obs_outcome)
  doubled_data_theta = rbind(data1, data2)
  
  w <- doubled_data_theta$w
  
  theta_update = tryCatch({glm(y ~ . -w -y, weights = w, data = doubled_data_theta,
                     family = "binomial"(link = "logit"))},
                     error = function(e){
                       glm(y ~ . -w -y,
                           #weights = w,
                           data = doubled_data_theta,
                           family = "binomial"(link = "logit"))
                     })

  theta_new <- unname(stats::coef(theta_update))
  
  # Save new parameters
  param_new = c(beta_new, gamma1_new, gamma2_new, theta_new)
  param_current = param_new
  return(param_new)

}

