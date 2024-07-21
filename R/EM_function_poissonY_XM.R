#' EM Algorithm Function for Estimation of the Misclassification Model
#' 
#' Function is for cases with \eqn{Y \sim Poisson} and with an interaction term
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
#' @importFrom stats coefficients binomial optim coef poisson
#'

EM_function_poissonY_XM <- function(param_current,
                                obs_mediator, obs_outcome,
                                X, Z, c_matrix,
                                sample_size, n_cat){

  # Create design matrix for true mediator model
  design_matrix = cbind(X, c_matrix)
  
  # Set up parameter indices
  gamma_index_1 = ncol(design_matrix) + 1
  gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
  
  n_param <- length(param_current)
  
  beta_current = matrix(param_current[1:ncol(design_matrix)], ncol = 1)
  gamma_current = matrix(c(param_current[gamma_index_1:gamma_index_2]),
                         ncol = n_cat, byrow = FALSE)
  theta_current = matrix(c(param_current[(gamma_index_2 + 1):n_param]),
                         ncol = 1)

  probabilities = pi_compute(beta_current, design_matrix, sample_size, n_cat)
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)
  
  # Compute likelihood value of Y based on x, m, c, theta, and sigma
  interaction_term_m0 <- X[,-1] * 0
  outcome_design_matrix_m0 <- cbind(cbind(X, cbind(rep(0, sample_size), c_matrix)),
                                    interaction_term_m0)
  model_y_m0 <- outcome_design_matrix_m0 %*% theta_current
  lambda_y_m0 = exp(model_y_m0)
  p_y_m0_term = ((lambda_y_m0 ^ obs_outcome) * exp(-lambda_y_m0)) / (factorial(obs_outcome))
  
  interaction_term_m1 <- X[,-1] * 1
  outcome_design_matrix_m1 <- cbind(cbind(X, cbind(rep(1, sample_size), c_matrix)),
                                    interaction_term_m1)
  model_y_m1 <- outcome_design_matrix_m1 %*% theta_current
  lambda_y_m1 = exp(model_y_m1)
  p_y_m1_term = ((lambda_y_m1 ^ obs_outcome) * exp(-lambda_y_m1)) / (factorial(obs_outcome))

  mstar_matrix = matrix(c(ifelse(obs_mediator == 1, 1, 0), 
                          ifelse(obs_mediator == 2, 1, 0)),
                        nrow = sample_size, byrow = FALSE)
  outcome_matrix = matrix(c(obs_outcome,
                            1 - obs_outcome),
                          nrow = sample_size, byrow = FALSE)
  weights = w_m_poissonY(mstar_matrix, outcome_matrix,
                        pistar_matrix = conditional_probabilities,
                        pi_matrix = probabilities,
                        p_y_m0_term, p_y_m1_term,
                        sample_size, n_cat)

  Mstar01 = mstar_matrix[,1]
  fit.gamma1 <- suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,1],
                           family = "binomial"(link = "logit")) )
  gamma1_new <- unname(coefficients(fit.gamma1))

  fit.gamma2 <- suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,2],
                           family = "binomial"(link = "logit")) )
  gamma2_new <- unname(coefficients(fit.gamma2))

  fit.beta <- suppressWarnings( stats::glm(weights[,1] ~ . + 0, as.data.frame(design_matrix),
                         family = stats::binomial()) )
  beta_new <- unname(coefficients(fit.beta))

  x_vector = X[,2]
  data1 = data.frame(x = x_vector, m = 0, c = c_matrix,
                     xm = x_vector * 0,
                     w = weights[,2], y = obs_outcome)
  data2 = data.frame(x = x_vector, m = 1, c = c_matrix, 
                     xm = x_vector * 1,
                     w = weights[,1], y = obs_outcome)
  doubled_data_theta = rbind(data1, data2)
  
  theta_update = glm(y ~ . -w -y, weights = doubled_data_theta$w,
                     data = doubled_data_theta,
                     family = "poisson"(link = "log"))

  theta_new <- unname(coef(theta_update))
  
  param_new = c(beta_new, gamma1_new, gamma2_new, theta_new)
  param_current = param_new
  return(param_new)

}

