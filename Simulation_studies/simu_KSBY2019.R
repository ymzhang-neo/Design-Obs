# Simulation: KSBY2019
# Kunzel, Sekhon, Bickel, Yu (2019)
# Metalearners for estimating heterogeneous treatment effects using machine learning
# Proceedings of the National Academy of Sciences of the United States of America

set.seed(1) 

# Notes -------------------------------------------------------------------
#
# [TODO] Implement the vine method of simulating correlation matrices.
#   Currently the correlation matrix is the identity matrix.

# 0 Covariates and setting up ---------------------------------------------

# # We do not consider separating the training and testing data, for now.
# para_N_train = 100000
# para_N_test = 100000
# para_N_total = para_N_train + para_N_test
# 
# dat_X_trainID = sample(1:para_N_total, para_N_train)
# dat_X_train = dat_X[dat_X_trainID]
# dat_X_test = dat_X[-dat_X_trainID]

# d: defined in each scenario

fn_getX = function(N, d, Sigma = NULL) {
  if (is.null(Sigma)) {
    Sigma = diag(d)
  }
  dat_X = t(t(chol(Sigma)) %*% matrix(rnorm(N * d), nrow = d))
  return(dat_X)
}

fn_getW = function(dat_X, fn_ps) {
  ps_prob = apply(dat_X, 1, fn_ps)
  dat_W = rbinom(n = nrow(dat_X), size = 1, prob = ps_prob)
  return(dat_W)
}


# 1 Unbalanced, simple CATE -----------------------------------------------

para_1_d = 20

fn_1_ps = function(x) {0.01}

para_1_mu_w0_beta = runif(para_1_d, min = -5, max = 5)

fn_1_mu_w0 = function(x) {
  sum(para_1_mu_w0_beta * x) + 5 * (x[1] > 0.5)
}

fn_1_mu_w1 = function(x) {
  sum(para_1_mu_w0_beta * x) + 5 * (x[1] > 0.5) + 8 * (x[2] > 0.1)
}


# 2 Balanced, no confounding, complex linear CATE -------------------------

para_2_d = 20

fn_2_ps = function(x) {0.5}

para_2_mu_w0_beta = runif(para_2_d, min = 1, max = 30)

para_2_mu_w1_beta = runif(para_2_d, min = 1, max = 30)

fn_2_mu_w0 = function(x) {sum(para_2_mu_w0_beta * x)}

fn_2_mu_w1 = function(x) {sum(para_2_mu_w1_beta * x)}


# 3 Balanced, no confounding, complex linear CATE -------------------------
# Ref: Wager and Athey (2018)

para_3_d = 20

fn_3_ps = function(x) {0.5}

fn_varsigma = function(x) {2 / (1 + exp(-12 * (x - 0.5)))}

fn_3_mu_w0 = function(x) {fn_varsigma(x[1]) * fn_varsigma(x[2]) / 2}

fn_3_mu_w1 = function(x) {- fn_varsigma(x[1]) * fn_varsigma(x[2]) / 2}


# 4 Balanceed, no confounding, no TE, linear mean -------------------------

para_4_d = 5

fn_4_ps = function(x) {0.5}

para_4_mu_beta = runif(para_4_d, min = 1, max = 30)

fn_4_mu_w0 = function(x) {sum(para_4_mu_beta * x)} 

fn_4_mu_w1 = fn_4_mu_w0

# 5 Balanceed, no confounding, no TE, piecewise linear mean ---------------

para_5_d = 20

fn_5_ps = function(x) {0.5}

para_5_mu_beta = runif(para_5_d, min = -15, max = 15)

fn_5_mu_w0 = function(x) {
  if (x[20] < -0.4) {
    nonzero_label = c(rep(1, 5), rep(0, para_5_d - 5))  
  } else if (x[20] <= 0.4) {
    nonzero_label = c(rep(0, 5), rep(1, 5), rep(0, para_5_d - 10))
  } else if (x[20] > 0.4) {
    nonzero_label = c(rep(0, 10), rep(1, 5), rep(0, para_5_d - 15))
  }
  mu = sum(para_5_mu_beta * x * nonzero_label)
  return(mu)
}

fn_5_mu_w1 = fn_5_mu_w0


# 6 Confounding -----------------------------------------------------------

para_6_d = 20

fn_6_ps = function(x) {(1 + dbeta(x[1], shape1 = 2, shape2 = 4)) / 4}

fn_6_mu_w0 = function(x) {2 * x[1] - 1}

fn_6_mu_w1 = fn_6_mu_w0

# Get the data! -----------------------------------------------------------

para_N = 10000 # 20k, 30k, 50k, 100k, 200k, 300k

para_d = para_1_d # para_2_d, etc.

fn_ps = fn_1_ps # fn_2_ps, etc.

fn_mu_w0 = fn_1_mu_w0 # fn_2_mu_w0, etc.

fn_mu_w1 = fn_1_mu_w1

dat_X = fn_getX(N = para_N, d = para_d)
# dat_X = matrix(runif(para_N * para_6_d), ncol = para_6_d) # Case 6 only

dat_W = fn_getW(dat_X = dat_X, fn_ps = fn_ps)

dat_PO = cbind(apply(dat_X, 1, fn_mu_w0), apply(dat_X, 1, fn_mu_w1)) + 
  matrix(rnorm(para_N * 2), ncol = 2)

dat_Y = dat_PO[, 1] * (1 - dat_W) + dat_PO[, 2] * dat_W