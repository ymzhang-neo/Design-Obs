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
para_N_train = 100000
para_N_test = 100000
para_N_total = para_N_train + para_N_test

para_d = 10
para_Sigma = diag(para_d)
para_Mu = rep(0, para_d)

dat_X = t(
  para_Mu + t(chol(para_Sigma)) %*% 
  matrix(rnorm(para_N_total * para_d), nrow = para_d)
)

dat_X_trainID = sample(1:para_N_total, para_N_train)
dat_X_train = dat_X[dat_X_trainID]
dat_X_test = dat_X[-dat_X_trainID]