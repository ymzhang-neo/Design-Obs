# Simulation: SWYYC2015
# Sun, Wang, Yin, Yang, Chang (2015)
# Causal inference via sparse additive models with application to online advertising
# AAAI'15 / Proceedings of the Twenty-Ninth AAAI Conference on Artificial Intelligence
# Also used in LVKF2016, LF2017

# Notes -------------------------------------------------------------------
# - The simulation is simplified by having only the binary treatment factor.
#   The original paper considered a three-dimensional treatment.

set.seed(20190918)

para_N = 1000
para_p = 200

fn_f1 = function(x) {return( -2 * sin(2 * x) )}
fn_f2 = function(x) {return( x ** 2 - 1 / 3 )}
fn_f3 = function(x) {return( x - 0.5 )}
fn_f4 = function(x) {return( exp(-x) - exp(-1) - 1 )}
fn_f5 = function(x) {return( (x - 0.5) ** 2 + 2 )}
fn_f6 = function(x) {return( as.numeric(x > 0) )}
fn_f7 = function(x) {return( exp(-x) )}
fn_f8 = function(x) {return( cos(x) )}
fn_f9 = function(x) {return( x ** 2 )}
fn_f10 = function(x) {return( x )}

list_fn = list(fn_f1, fn_f2, fn_f3, fn_f4, fn_f5, fn_f6, fn_f7, fn_f8, fn_f9, fn_f10)


dat_X = matrix(rnorm(para_N * para_p), nrow = para_N)

dat_W_latent = sapply(dat_X[, 1], list_fn[[1]]) + 
               sapply(dat_X[, 2], list_fn[[2]]) + 
               sapply(dat_X[, 3], list_fn[[3]]) + 
               sapply(dat_X[, 4], list_fn[[4]]) + 
               sapply(dat_X[, 5], list_fn[[5]])
dat_W = as.numeric(dat_W_latent > 0)

dat_Y_mean = dat_W + 
             sapply(dat_X[, 1], list_fn[[6]]) + 
             sapply(dat_X[, 2], list_fn[[7]]) + 
             sapply(dat_X[, 3], list_fn[[8]]) + 
             sapply(dat_X[, 4], list_fn[[9]]) + 
             sapply(dat_X[, 5], list_fn[[10]])
dat_Y = rnorm(para_N) + dat_Y_mean
