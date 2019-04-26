# Simulation: LD2004
# Lunceford, Davidian (2004)
# Stratification and weighting via the propensity score in estimation of causal treatment effects: a comparative study
# Statistics in Medicine

set.seed(1)

para_N = 1000 # 5000

para_mu_X1V1X2V2_1 = c(1, 1, -1, -1)
para_mu_X1V1X2V2_0 = c(-1, -1, 1, 1)
para_cov_X1V1X2V2 = matrix(
  c(  1,  .5, -.5, -.5, 
     .5,   1, -.5, -.5,
    -.5, -.5,   1,  .5,
    -.5, -.5,  .5,   1),
  nrow = 4, ncol = 4
)

para_beta_W = c(0, 0.6, -0.6, 0.6) # str
# para_beta_W = c(0, 0.3, -0.3, 0.3) # mod

para_beta_Y_W = 2
para_beta_Y_X = c(0, -1, 1, -1)
para_beta_Y_V = c(-1, 1, 1) # str
# para_beta_Y_V = c(-0.5, 0.5, 0.5) # mod
# para_beta_Y_V = c(0, 0, 0) # no
para_epsilon_sd = 1


dat_X3 = rbinom(para_N, size = 1, prob = 0.2)
dat_V3 = rbinom(para_N, size = 1, prob = 0.75*dat_X3 + 0.25*(1-dat_X3))

dat_X1V2X2V2 = t(t(chol(para_cov_X1V1X2V2)) %*% matrix(rnorm(para_N * 4), nrow = 4)) + 
  diag(dat_X3) %*% matrix(para_mu_X1V1X2V2_1, byrow = TRUE, nrow = para_N, ncol = 4) + 
  diag(abs(dat_X3 - 1)) %*% matrix(para_mu_X1V1X2V2_0, byrow = TRUE, nrow = para_N, ncol = 4)
  # diag(dat_X3 - 1) %*% matrix(para_mu_X1V1X2V2_1, byrow = TRUE, nrow = para_N, ncol = 4)

dat_X = cbind(dat_X1V2X2V2, dat_X3, dat_V3)
colnames(dat_X) = c('X1', 'V1', 'X2', 'V2', 'X3', 'V3')


dat_W = rbinom(para_N, size = 1, 
               prob = 1/(1 + exp(- (dat_X[, c('X1', 'X2', 'X3')] %*% para_beta_W[-1] + para_beta_W[1]))))


# Simulate outcomes -------------------------------------------------------

dat_Y = para_beta_Y_W * dat_W + 
  para_beta_Y_X[1] + dat_X[, c('X1', 'X2', 'X3')] %*% para_beta_Y_X[-1] + 
  dat_X[, c('V1', 'V2', 'V3')] %*% para_beta_Y_V + 
  rnorm(para_N, sd = para_epsilon_sd)