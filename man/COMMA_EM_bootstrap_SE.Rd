% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMMA_EM_bootstrap_SE.R
\name{COMMA_EM_bootstrap_SE}
\alias{COMMA_EM_bootstrap_SE}
\title{Estimate Bootstrap Standard Errors using EM}
\usage{
COMMA_EM_bootstrap_SE(
  parameter_estimates,
  sigma_estimate = 1,
  n_bootstrap,
  n_parallel,
  outcome_distribution,
  interaction_indicator,
  x_matrix,
  z_matrix,
  c_matrix,
  tolerance = 1e-07,
  max_em_iterations = 1500,
  em_method = "squarem",
  random_seed = NULL
)
}
\arguments{
\item{parameter_estimates}{A column matrix of \eqn{\beta}, \eqn{\gamma},
and \eqn{\theta} parameter values obtained from a COMMA analysis function.
Parameter estimates should be supplied in the following order: 1) \eqn{\beta}
(intercept, slope), 2) \eqn{\gamma} (intercept and slope from the M = 1
mechanism, intercept and slope from the M = 2 mechanism), and 3) \eqn{\theta}
(intercept, slope, coefficient for \code{x}, slope coefficient for \code{m},
slope coefficient for \code{c}, and, optionally, slope coefficient for
\code{xm} if using).}

\item{sigma_estimate}{A numeric value specifying the estimated
standard deviation. This value is only required if \code{outcome_distribution}
is \code{"Normal"}. Default is 1. For non-Normal outcome distributions, the
value should be \code{NULL}.}

\item{n_bootstrap}{A numeric value specifying the number of bootstrap samples
to draw.}

\item{n_parallel}{A numeric value specifying the number of parallel cores to
run the computation on.}

\item{outcome_distribution}{A character string specifying the distribution of
the outcome variable. Options are \code{"Bernoulli"}, \code{"Normal"}, or
\code{"Poisson"}.}

\item{interaction_indicator}{A logical value indicating if an interaction between
\code{x} and \code{m} should be used to generate the outcome variable, \code{y}.}

\item{x_matrix}{A numeric matrix of predictors in the true mediator and outcome mechanisms.
\code{x_matrix} should not contain an intercept and no values should be \code{NA}.}

\item{z_matrix}{A numeric matrix of covariates in the observation mechanism.
\code{z_matrix} should not contain an intercept and no values should be \code{NA}.}

\item{c_matrix}{A numeric matrix of covariates in the true mediator and outcome mechanisms.
\code{c_matrix} should not contain an intercept and no values should be \code{NA}.}

\item{tolerance}{A numeric value specifying when to stop estimation, based on
the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.}

\item{max_em_iterations}{A numeric value specifying when to stop estimation, based on
the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.}

\item{em_method}{A character string specifying which EM algorithm will be applied.
Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
recommended option is \code{"squarem"}.}

\item{random_seed}{A numeric value specifying the random seed to set for bootstrap
sampling. Default is \code{NULL}.}
}
\value{
\code{COMMA_EM_bootstrap_SE} returns a list with two elements: 1)
\code{bootstrap_df} and 2) \code{bootstrap_SE}. \code{bootstrap_df} is a data
frame containing \code{COMMA_EM} output for each bootstrap sample. \code{bootstrap_SE}
is a data frame containing bootstrap standard error estimates for each parameter.
}
\description{
Estimate Bootstrap Standard Errors using EM
}
\examples{
\donttest{
set.seed(20240709)
sample_size <- 2000

n_cat <- 2 # Number of categories in the binary mediator

# Data generation settings
x_mu <- 0
x_sigma <- 1
z_shape <- 1
c_shape <- 1

# True parameter values (gamma terms set the misclassification rate)
true_beta <- matrix(c(1, -2, .5), ncol = 1)
true_gamma <- matrix(c(1, 1, -.5, -1.5), nrow = 2, byrow = FALSE)
true_theta <- matrix(c(1, 1.5, -2, -.2), ncol = 1)

example_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
                           interaction_indicator = FALSE,
                           outcome_distribution = "Bernoulli",
                           true_beta, true_gamma, true_theta)
                           
beta_start <- matrix(rep(1, 3), ncol = 1)
gamma_start <- matrix(rep(1, 4), nrow = 2, ncol = 2)
theta_start <- matrix(rep(1, 4), ncol = 1)

Mstar = example_data[["obs_mediator"]]
outcome = example_data[["outcome"]]
x_matrix = example_data[["x"]]
z_matrix = example_data[["z"]]
c_matrix = example_data[["c"]]
                           
EM_results <- COMMA_EM(Mstar, outcome, "Bernoulli", FALSE,
                       x_matrix, z_matrix, c_matrix,
                       beta_start, gamma_start, theta_start)

EM_results

EM_SEs <- COMMA_EM_bootstrap_SE(EM_results$Estimates, sigma_estimate = NULL,
                                n_bootstrap = 3,
                                n_parallel = 1,
                                outcome_distribution = "Bernoulli",
                                interaction_indicator = FALSE,
                                x_matrix, z_matrix, c_matrix,
                                random_seed = 1)                
EM_SEs$bootstrap_SE
}
}
