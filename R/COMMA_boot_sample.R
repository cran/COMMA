#' Generate Bootstrap Samples for Estimating Standard Errors
#'
#' @param parameter_estimates A column matrix of \eqn{\beta}, \eqn{\gamma},
#' and \eqn{\theta} parameter values obtained from a COMMA analysis function.
#' Parameter estimates should be supplied in the following order: 1) \eqn{\beta}
#' (intercept, slope), 2) \eqn{\gamma} (intercept and slope from the M = 1 
#' mechanism, intercept and slope from the M = 2 mechanism), and 3) \eqn{\theta}
#' (intercept, slope, coefficient for \code{x}, slope coefficient for \code{m},
#' slope coefficient for \code{c}, and, optionally, slope coefficient for
#' \code{xm} if using).
#' @param sigma_estimate A numeric value specifying the estimated 
#'   standard deviation. This value is only required if \code{outcome_distribution}
#'   is \code{"Normal"}. Default is 1. For non-Normal outcome distributions, the 
#'   value should be \code{NULL}.
#' @param outcome_distribution A character string specifying the distribution of 
#'   the outcome variable. Options are \code{"Bernoulli"}, \code{"Normal"}, or
#'   \code{"Poisson"}.
#' @param interaction_indicator A logical value indicating if an interaction between
#'   \code{x} and \code{m} should be used to generate the outcome variable, \code{y}.
#' @param x_matrix A numeric matrix of predictors in the true mediator and outcome mechanisms.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param c_matrix A numeric matrix of covariates in the true mediator and outcome mechanisms.
#'   \code{c_matrix} should not contain an intercept and no values should be \code{NA}.
#'
#' @return \code{COMMA_boot_sample} returns a list with the bootstrap sample data:
#'   \item{obs_mediator}{A vector of observed mediator values.}
#'   \item{true_mediator}{A vector of true mediator values.}
#'   \item{outcome}{A vector of outcome values.}
#'   \item{x_matrix}{A matrix of predictor values in the true mediator
#'   mechanism. Identical to that supplied by the user.}
#'   \item{z_matrix}{A matrix of predictor values in the observed mediator
#'   mechanism. Identical to that supplied by the user.}
#'   \item{c_matrix}{A matrix of covariates. Identical to that supplied by the user.}
#'   
#' @include pi_compute.R
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom rpois
#'
COMMA_boot_sample <- function(parameter_estimates,
                              sigma_estimate = 1,
                              outcome_distribution,
                              interaction_indicator,
                              # Predictor matrices
                              x_matrix, z_matrix, c_matrix){
  
  n_cat = 2 # Number of categories in mediator
  sample_size = length(c(x_matrix)) # Sample size
  
  # Create design matrices
  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)
  
  x_mat <- as.matrix(x_matrix)
  c_mat <- as.matrix(c_matrix)
  
  # Create matrix of true mediation model predictors
  mediation_model_predictors <- cbind(x_matrix, c_matrix)
  
  beta_param <- parameter_estimates[1:(ncol(mediation_model_predictors) + 1)]
  gamma_param <- matrix(parameter_estimates[(ncol(mediation_model_predictors) + 2):(
    ncol(mediation_model_predictors) + 2 + ((ncol(Z) * 2)) - 1)],
    ncol = 2, byrow = FALSE)
  theta_param <- parameter_estimates[(2 + ((ncol(mediation_model_predictors) + (ncol(Z) * 2)))):length(parameter_estimates)]
  
  predictor_matrix <- cbind(1, mediation_model_predictors)
  # Generate probabilities for the true mediator value
  pi_matrix <- pi_compute(beta_param, predictor_matrix, sample_size, n_cat)
  
  # Generate true mediator variable based on probabilities in pi_matrix
  true_M <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_M[i] = which(rmultinom(1, 1, pi_matrix[i,]) == 1)
  }
  
  # Generate probabilities for observed mediator conditional on true mediator
  pistar_matrix <- pistar_compute(gamma_param, Z, sample_size, n_cat)
  
  # Generate observed mediator variable based on conditional probabilities in pistar_matrix
  obs_Mstar <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_j = true_M[i]
    obs_Mstar[i] = which(rmultinom(1, 1,
                                   pistar_matrix[c(i,sample_size + i),
                                                 true_j]) == 1)
  }
  
  # Create a matrix of observed mediator variables using dummy coding
  obs_Mstar_reps <- matrix(rep(obs_Mstar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Mstar_matrix <- 1 * (obs_Mstar_reps == category_matrix)
  
  # Indicator variable for true mediator.
  m_indicator <- ifelse(true_M == 1, 1, 0)
  interaction_term <- m_indicator * x_matrix
  
  if(interaction_indicator == FALSE & outcome_distribution == "Bernoulli"){
    
    # Generate probabilities for the outcome
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x_matrix, m_indicator, c_matrix),
                                    nrow = sample_size, byrow = FALSE)
    exp_result <- exp(outcome_design_matrix %*% theta_param)
    
    exp_denominator <- 1 + exp_result
    
    classification_prob_matrix <- matrix(c(exp_result / exp_denominator,
                                           1 / exp_denominator),
                                         ncol = 2, byrow = FALSE)
    
    # Generate binary outcome based on probabilities in pitilde_matrix
    outcome_12 <- rep(NA, sample_size)
    for(i in 1:sample_size){
      outcome_12[i] = which(rmultinom(1, 1,
                                      classification_prob_matrix[i,]) == 1)
    }
    
    # Convert outcome to 0/1 instead of 2/1 variable coding
    outcome = ifelse(outcome_12 == 1, 1, 0)
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Bernoulli"){
    
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x_matrix, m_indicator, c_matrix,
                                      interaction_term),
                                    nrow = sample_size, byrow = FALSE)
    exp_result <- exp(outcome_design_matrix %*% theta_param)
    
    exp_denominator <- 1 + exp_result
    
    classification_prob_matrix <- matrix(c(exp_result / exp_denominator,
                                           1 / exp_denominator),
                                         ncol = 2, byrow = FALSE)
    
    # Generate binary outcome based on probabilities in pitilde_matrix
    outcome_12 <- rep(NA, sample_size)
    for(i in 1:sample_size){
      outcome_12[i] = which(rmultinom(1, 1,
                                      classification_prob_matrix[i,]) == 1)
    }
    
    # Convert outcome to 0/1 instead of 2/1 variable coding
    outcome = ifelse(outcome_12 == 1, 1, 0)
    
  } else if(interaction_indicator == FALSE & outcome_distribution == "Normal"){
    
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x_matrix, m_indicator, c_matrix),
                                    nrow = sample_size, byrow = FALSE)
    
    # Generate mean and Normal errors
    additive_term <- outcome_design_matrix %*% theta_param
    additive_term_error <- rnorm(sample_size, sigma_estimate) 
    
    # Return value of Y
    outcome <- additive_term + additive_term_error
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Normal"){
    
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x_matrix, m_indicator, c_matrix,
                                      interaction_term),
                                    nrow = sample_size, byrow = FALSE)
    
    # Generate mean and Normal errors
    additive_term <- outcome_design_matrix %*% theta_param
    additive_term_error <- rnorm(sample_size, sigma_estimate) 
    
    # Return value of Y
    outcome <- additive_term + additive_term_error
    
  } else if(interaction_indicator == FALSE & outcome_distribution == "Poisson"){
    
    # Generate probabilities for the outcome
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x_matrix, m_indicator, c_matrix),
                                    nrow = sample_size, byrow = FALSE)
    
    # Generate mean
    exp_result <- exp(outcome_design_matrix %*% theta_param)
    
    # Return value of Y
    outcome <- rpois(sample_size, exp_result)
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Poisson"){
    
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x_matrix, m_indicator, c_matrix,
                                      interaction_term),
                                    nrow = sample_size, byrow = FALSE)
    # Generate mean
    exp_result <- exp(outcome_design_matrix %*% theta_param)
    
    # Return value of Y
    outcome <- rpois(sample_size, exp_result)
    
    
  } else {
    outcome = "Undefined"
  }
  
  # Organize data for output
  data_output <- list(obs_mediator = obs_Mstar,
                      outcome = outcome,
                      true_mediator = true_M,
                      x_matrix = x_matrix,
                      z_matrix = z_matrix,
                      c_matrix = c_matrix)
  
  return(data_output)
  
}
