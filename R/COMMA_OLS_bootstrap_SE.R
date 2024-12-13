#' Estimate Bootstrap Standard Errors using OLS
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
#'   standard deviation. Default is 1. 
#' @param n_bootstrap A numeric value specifying the number of bootstrap samples
#' to draw. 
#' @param n_parallel A numeric value specifying the number of parallel cores to
#' run the computation on.
#' @param x_matrix A numeric matrix of predictors in the true mediator and outcome mechanisms.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param c_matrix A numeric matrix of covariates in the true mediator and outcome mechanisms.
#'   \code{c_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param em_method A character string specifying which EM algorithm will be applied.
#'   Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
#'   recommended option is \code{"squarem"}.
#' @param random_seed A numeric value specifying the random seed to set for bootstrap
#'   sampling. Default is \code{NULL}.
#'
#' @return \code{COMMA_OLS_bootstrap_SE} returns a list with two elements: 1)
#' \code{bootstrap_df} and 2) \code{bootstrap_SE}. \code{bootstrap_df} is a data
#' frame containing \code{COMMA_OLS} output for each bootstrap sample. \code{bootstrap_SE}
#' is a data frame containing bootstrap standard error estimates for each parameter. 
#' 
#' @export
#' 
#' @include COMMA_OLS.R
#' @include COMMA_boot_sample.R
#' 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% 
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr %>% group_by ungroup summarise
#' @importFrom stats sd
#'
#' @examples \donttest{
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
#' true_theta <- matrix(c(1, 1.5, -2, 2), ncol = 1)
#' 
#' example_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
#'                            interaction_indicator = FALSE,
#'                            outcome_distribution = "Normal",
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
#' OLS_results <- COMMA_OLS(Mstar, outcome,
#'                          x_matrix, z_matrix, c_matrix,
#'                          beta_start, gamma_start, theta_start)
#'                          
#' OLS_results
#' 
#' OLS_SEs <- COMMA_OLS_bootstrap_SE(OLS_results$Estimates, sigma_estimate = 1,
#'                                   n_bootstrap = 3,
#'                                   n_parallel = 1,
#'                                   x_matrix, z_matrix, c_matrix,
#'                                   random_seed = 1)
#'                                   
#' OLS_SEs$bootstrap_SE
#' }
COMMA_OLS_bootstrap_SE <- function(parameter_estimates,
                                   sigma_estimate = 1,
                                   n_bootstrap,
                                   n_parallel,
                                   # Predictor matrices
                                   x_matrix, z_matrix, c_matrix,
                                   tolerance = 1e-7,
                                   max_em_iterations = 1500,
                                   em_method = "squarem",
                                   random_seed = NULL){
  
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
  
  beta_start <- parameter_estimates[1:(ncol(mediation_model_predictors) + 1)]
  gamma_start <- matrix(parameter_estimates[(ncol(mediation_model_predictors) + 2):(
    ncol(mediation_model_predictors) + 2 + ((ncol(Z) * 2)) - 1)],
    ncol = 2, byrow = FALSE)
  theta_start <- parameter_estimates[(2 + ((ncol(mediation_model_predictors) + (ncol(Z) * 2)))):length(parameter_estimates)]
  
  i = NULL
  Parameter = NULL
  Estimates = NULL
  
  cluster <- parallel::makeCluster(n_parallel) 
  doParallel::registerDoParallel(cluster)
  
  bootstrap_df <- foreach(i = 1:n_bootstrap,
          .combine = rbind) %dopar% {
            
    if(is.null(random_seed)){
      random_seed_i <- NULL
    } else {
      random_seed_i <- random_seed + i
    }
    
    set.seed(random_seed_i)
    
    boot_sample_i <- COMMA_boot_sample(parameter_estimates,
                                       sigma_estimate,
                                       outcome_distribution = "Normal",
                                       interaction_indicator = FALSE,
                                       # Predictor matrices
                                       x_matrix, z_matrix, c_matrix)
    
    bootstrap_param <- COMMA_OLS(boot_sample_i[["obs_mediator"]],
                                 boot_sample_i[["outcome"]],
                                 x_matrix, z_matrix, c_matrix,
                                 beta_start, gamma_start, theta_start,
                                 tolerance, max_em_iterations, em_method)
    
    bootstrap_param$bootstrap_iteration <- i
    
    bootstrap_param
    
          }
  
  stopCluster(cluster)
  
  bootstrap_SE <- bootstrap_df %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(Mean = mean(Estimates, na.rm = TRUE),
                     SE = stats::sd(Estimates, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  bootstrap_results <- list(bootstrap_df = bootstrap_df,
                            bootstrap_SE = bootstrap_SE)
  
  
}