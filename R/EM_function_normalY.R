#' EM Algorithm Function for Estimation of the Misclassification Model
#' 
#' Function is for cases with \eqn{Y \sim Normal} and with no interaction term
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
#' @importFrom stats rnorm rgamma rmultinom coefficients binomial optim
#'
EM_function_normalY <- function(param_current,
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
  theta_current = matrix(c(param_current[(gamma_index_2 + 1):(n_param - 1)]),
                         ncol = 1)
  sigma_current = param_current[n_param]
  
  # Compute probability of each latent mediator value
  probabilities = pi_compute(beta_current, design_matrix, sample_size, n_cat)
  
  # Compute probability of observed mediator, given latent mediator
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)
  
  # Compute likelihood value of Y based on x, m, c, theta, and sigma
  outcome_design_matrix_m0 <- cbind(X, cbind(rep(0, sample_size), c_matrix))
  model_y_m0 <- outcome_design_matrix_m0 %*% theta_current
  residual_term_m0 = obs_outcome - model_y_m0
  term1_m0 = 1 / sqrt(2 * pi * c(sigma_current ^ 2))
  exp_term_m0 = exp(-1 * residual_term_m0^2 * (1 / c(2 * sigma_current^2)))
  p_yi_m0 = term1_m0 * exp_term_m0
  
  outcome_design_matrix_m1 <- cbind(X, cbind(rep(1, sample_size), c_matrix))
  model_y_m1 <- outcome_design_matrix_m1 %*% theta_current
  residual_term_m1 = obs_outcome - model_y_m1
  term1_m1 = 1 / sqrt(2 * pi * c(sigma_current ^ 2))
  exp_term_m1 = exp(-1 * residual_term_m1^2 * (1 / c(2 * sigma_current^2)))
  p_yi_m1 = term1_m1 * exp_term_m1

  # Create a matrix of observed mediator variables using dummy coding
  mstar_matrix = matrix(c(ifelse(obs_mediator == 1, 1, 0), 
                          ifelse(obs_mediator == 2, 1, 0)),
                        nrow = sample_size, byrow = FALSE)
  
  # Compute E-Step weights
  weights = w_m_normalY(mstar_matrix,
                        pistar_matrix = conditional_probabilities,
                        pi_matrix = probabilities,
                        p_yi_m0, p_yi_m1,
                        sample_size, n_cat)

  # Estimate gamma parameters using weighted logistic regression
  ## Weights from E-Step (split by value of latent mediator, m)
  ## Outcome is the observed mediator
  Mstar01 = mstar_matrix[,1]
  fit.gamma1 <- suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,1],
                           family = "binomial"(link = "logit")) )
  gamma1_new <- unname(coefficients(fit.gamma1))

  fit.gamma2 <- suppressWarnings( stats::glm(Mstar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,2],
                           family = "binomial"(link = "logit")) )
  gamma2_new <- unname(coefficients(fit.gamma2))

  # Estimate beta parameters using logistic regression
  ## Outcome is the E-Step weight
  fit.beta <- suppressWarnings( stats::glm(weights[,1] ~ . + 0, as.data.frame(design_matrix),
                         family = stats::binomial()) )
  beta_new <- unname(coefficients(fit.beta))

  if(ncol(matrix(c_matrix, nrow = sample_size, byrow = FALSE)) == 1){
    
    # Solve for theta parameters using a system of equations
    a_row1 <- c(1, mean(X[,2]), mean(c_matrix), mean(weights[,1]))
    a_row2 <- c(sum(X[,2]) / sum(X[,2]^2), 1, sum(X[,2] * c_matrix) / sum(X[,2]^2),
                sum(X[,2] * weights[,1]) / sum(X[,2]^2))
    a_row3 <- c(sum(c_matrix) / sum(c_matrix^2), sum(c_matrix * X[,2]) / sum(c_matrix^2),
                1, sum(c_matrix * weights[,1]) / sum(c_matrix^2))
    a_row4 <- c(1, sum(X[,2] * weights[,1]) / sum(weights[,1]),
                sum(c_matrix * weights[,1]) / sum(weights[,1]), 1)
    
    A = matrix(c(a_row1, a_row2, a_row3, a_row4), byrow = TRUE, nrow = 4)
    B = matrix(c(mean(obs_outcome), sum(X[,2] * obs_outcome) / sum(X[,2]^2),
                 sum(c_matrix * obs_outcome) / sum(c_matrix^2),
                 sum(weights[,1] * obs_outcome) / sum(weights[,1])),
               ncol = 1)
    
    theta_update <- solve(A, B)
    
    # Compute sigma estimate
    sigma_update <- (1 / sample_size) * sum(
      ((obs_outcome - theta_update[1] - theta_update[2] * X[,2]
        - theta_update[3] * c_matrix) ^ 2)
      - 2 * theta_update[4] * weights[,1] * (obs_outcome - theta_update[1] - theta_update[2] * X[,2]
                                             - theta_update[3] * c_matrix)
      + (theta_update[4] ^ 2 * weights[,1])
    )
    
    # Reorder theta estimates
    theta_new <- theta_update[c(1,2,4,3)]
    
    sigma_new <- sigma_update
    
  } else {
    
    theta_sigma_new_optim <- optim(par = c(theta_current, sigma_current),
                                   fn = theta_optim,
                                   m = weights[,1],
                                   x = X[,2], c_matrix = c_matrix,
                                   outcome = obs_outcome,
                                   sample_size = sample_size,
                                   n_cat = 2, method = "BFGS")
    
    theta_new <- theta_sigma_new_optim$par[-length(theta_sigma_new_optim$par)]
    sigma_new <- theta_sigma_new_optim$par[length(theta_sigma_new_optim$par)]
    
  }
  
  
  # Save new parameters
  param_new = c(beta_new, gamma1_new, gamma2_new, theta_new, sigma_new)
  param_new
  param_current = param_new
  return(param_new)

}




