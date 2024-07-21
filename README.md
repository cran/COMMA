# COMMA 

<!-- badges: start -->
![my workflow](https://github.com/kimberlywebb/COMMA/actions/workflows/r.yml/badge.svg)
<!-- badges: end -->

![ ](https://github.com/kimberlywebb/COMMA/blob/0c1c2c02f672d88543f2fc31a9813c242e24bfc6/COMMA_hex_sticker.png?raw=true)

**COMMA:** **CO**rrecting **M**isclassified **M**ediation **A**nalysis

Overview
--------------------------------------------------

**COMMA** provides a set of functions for the analysis of mediation models with binary mediator misclassification. 

The two main parts are:

- Classification probability calculations
- Parameter estimation 


Classification probability calculations
--------------------------------------------------
The package allows users to compute the probability of the latent true mediators and the conditional probability of observing a mediator value given the latent true mediator, based on parameters estimated from the `COMMA_EM`, `COMMA_PVW`, and `COMMA_OLS` functions.


Parameter estimation 
--------------------------------------------------
Jointly estimate parameters from the true mediator, observed mediator, and outcome mechanisms, respectively, in a mediation analysis with a binary misclassified mediator variable. Three methods are provided for parameter estimation: an ordinary least squares correction procedure, a predictive value weighting approach, and an expectation-maximization algorithm.

Installation
--------------------------------------------------

``` r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("kimberlywebb/COMMA")
```
