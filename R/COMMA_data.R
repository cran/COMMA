#' Generate Data to use in COMMA Functions
#'
#' @param sample_size An integer specifying the sample size of the generated data set.
#' @param x_mu A numeric value specifying the mean of \code{x} predictors
#' generated from a Normal distribution.
#' @param x_sigma A positive numeric value specifying the standard deviation of
#'   \code{x} predictors generated from a Normal distribution.
#' @param z_shape A positive numeric value specifying the shape parameter of
#'   \code{z} predictors generated from a Gamma distribution.
#' @param c_shape A positive numeric value specifying the shape parameter of
#'   \code{c} covariates generated from a Gamma distribution.
#' @param interaction_indicator A logical value indicating if an interaction between
#'   \code{x} and \code{m} should be used to generate the outcome variable, \code{y}.
#' @param outcome_distribution A character string specifying the distribution of 
#'   the outcome variable. Options are \code{"Bernoulli"}, \code{"Normal"}, or
#'   \code{"Poisson"}.
#' @param true_beta A column matrix of \eqn{\beta} parameter values (intercept, slope)
#'   to generate data under in the true mediator mechanism.
#' @param true_gamma A numeric matrix of \eqn{\gamma} parameters
#'   to generate data in the observed mediator mechanisms.
#'   In matrix form, the \code{gamma} matrix rows correspond to intercept (row 1)
#'   and slope (row 2) terms. The gamma parameter matrix columns correspond to the true mediator categories
#'   \eqn{M \in \{1, 2\}}.
#' @param true_theta A column matrix of \eqn{\theta} parameter values (intercept, slope
#'   coefficient for \code{x}, slope coefficient for \code{m}, slope coefficient for \code{c}),
#'   and, optionally, slope coefficient for \code{xm} if using) to generate data
#'   in the outcome mechanism.
#'
#' @return \code{COMMA_data} returns a list of generated data elements:
#'   \item{obs_mediator}{A vector of observed mediator values.}
#'   \item{true_mediator}{A vector of true mediator values.}
#'   \item{outcome}{A vector of outcome values.}
#'   \item{x}{A vector of generated predictor values in the true mediator
#'   mechanism, from the Normal distribution.}
#'   \item{z}{A vector of generated predictor values in the observed mediator
#'   mechanism from the Gamma distribution.}
#'   \item{c}{A vector of generated covariates.}
#'   \item{x_design_matrix}{The design matrix for the \code{x} predictor.}
#'   \item{z_design_matrix}{The design matrix for the \code{z} predictor.}
#'   \item{c_design_matrix}{The design matrix for the \code{c} predictor.}
#'
#' @export
#' 
#' @include pi_compute.R
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom rpois
#'
#' @examples 
#' set.seed(20240709)
#' sample_size <- 10000
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
#' head(example_data$obs_mediator)
#' head(example_data$true_mediator) 
#' 
COMMA_data <- function(sample_size,
                       x_mu, # Mean of X ~ Normal
                       x_sigma, # SD of X ~ Normal
                       z_shape, # Shape for Z ~ Gamma
                       c_shape, # Shape for C ~ Gamma
                       interaction_indicator, # Indicator for interaction term
                       outcome_distribution, # Distribution of outcome, Y
                       true_beta, true_gamma, true_theta){
  
  n_cat <- 2 # Number of categories in mediator
  
  # Generate X
  x <- rnorm(sample_size, x_mu, x_sigma)
  x_matrix <- matrix(c(rep(1, sample_size),
                       x),
                     nrow = sample_size, byrow = FALSE)

  # Generate Z
  z <- rgamma(sample_size, z_shape)
  z_matrix <- matrix(c(rep(1, sample_size),
                       z),
                     nrow = sample_size, byrow = FALSE)

  # Generate C
  c <- rgamma(sample_size, c_shape)
  c_matrix <- matrix(c(rep(1, sample_size),
                       c),
                     nrow = sample_size, byrow = FALSE)
  
  # Create matrix of predictors for the true mediator
  predictor_matrix <- cbind(x_matrix, c_matrix[,2])
  
  # Generate probabilities for the true mediator value
  pi_matrix <- pi_compute(true_beta, predictor_matrix, sample_size, n_cat)

  # Generate true mediator variable based on probabilities in pi_matrix
  true_M <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_M[i] = which(rmultinom(1, 1, pi_matrix[i,]) == 1)
  }
  
  # Generate probabilities for observed mediator conditional on true mediator
  pistar_matrix <- pistar_compute(true_gamma, z_matrix, sample_size, n_cat)
  
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
  interaction_term <- m_indicator * x

  if(interaction_indicator == FALSE & outcome_distribution == "Bernoulli"){
    
    # Generate probabilities for the outcome
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x, m_indicator, c),
                                    ncol = 4, byrow = FALSE)
    exp_result <- exp(outcome_design_matrix %*% true_theta)
    
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
                                      x, m_indicator, c,
                                      interaction_term),
                                    ncol = 5, byrow = FALSE)
    exp_result <- exp(outcome_design_matrix %*% true_theta)
    
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
                                      x, m_indicator, c),
                                    ncol = 4, byrow = FALSE)
    
    # Generate mean and Normal errors
    additive_term <- outcome_design_matrix %*% true_theta
    additive_term_error <- rnorm(sample_size) # Errors generated with SD, sigma = 1
    
    # Return value of Y
    outcome <- additive_term + additive_term_error
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Normal"){
    
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x, m_indicator, c,
                                      interaction_term),
                                    ncol = 5, byrow = FALSE)
    
    # Generate mean and Normal errors
    additive_term <- outcome_design_matrix %*% true_theta
    additive_term_error <- rnorm(sample_size) # Errors generated with SD, sigma = 1
    
    # Return value of Y
    outcome <- additive_term + additive_term_error
    
  } else if(interaction_indicator == FALSE & outcome_distribution == "Poisson"){
    
    # Generate probabilities for the outcome
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x, m_indicator, c),
                                    ncol = 4, byrow = FALSE)
    
    # Generate mean
    exp_result <- exp(outcome_design_matrix %*% true_theta)
    
    # Return value of Y
    outcome <- rpois(sample_size, exp_result)
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Poisson"){
    
    outcome_design_matrix <- matrix(c(rep(1, sample_size),
                                      x, m_indicator, c,
                                      interaction_term),
                                    ncol = 5, byrow = FALSE)
    # Generate mean
    exp_result <- exp(outcome_design_matrix %*% true_theta)
    
    # Return value of Y
    outcome <- rpois(sample_size, exp_result)
    
    
  } else {
    outcome = "Undefined"
  }
  
  
  # Organize data for output
  data_output <- list(obs_mediator = obs_Mstar,
                      outcome = outcome,
                      true_mediator = true_M,
                      x = x,
                      z = z,
                      c = c,
                      x_design_matrix = x_matrix,
                      z_design_matrix = z_matrix,
                      c_design_matrix = c_matrix)

  return(data_output)

}

