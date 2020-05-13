# Simulation: LD2004
# Lunceford, Davidian (2004)
# Stratification and weighting via the propensity score in estimation of causal treatment effects: a comparative study
# Statistics in Medicine

# Parameters --------------------------------------------------------------

n_sample = 1000 # 5000

para_X_fn = function(n_sample) {
  # input:
  #   n_sample: integer, sample size
  
  para_X_X3_p = 0.2 
  
  para_X_X1V1X2V2_mu_1 = c(1, 1, -1, -1)
  para_X_X1V1X2V2_mu_0 = c(-1, -1, 1, 1)
  para_X_X1V1X2V2_cov = matrix(
    c(  1,  .5, -.5, -.5, 
        .5,   1, -.5, -.5,
        -.5, -.5,   1,  .5,
        -.5, -.5,  .5,   1),
    nrow = 4, ncol = 4
  )
  
  dat_X_X3 = rbinom(n_sample, size = 1, prob = para_X_X3_p)
  
  dat_X_V3 = rbinom(n_sample, size = 1, 
                    prob = 0.75 * dat_X_X3 + 0.25 * (1 - dat_X_X3))
  
  dat_X_X1V1X2V2 = t(
    t(chol(para_X_X1V1X2V2_cov)) %*% 
      matrix(rnorm(n_sample * 4), nrow = 4)
  ) + 
    diag(dat_X_X3) %*% 
    matrix(para_X_X1V1X2V2_mu_1, byrow = TRUE, nrow = n_sample, ncol = 4) + 
    diag(1 - dat_X_X3) %*% 
    matrix(para_X_X1V1X2V2_mu_0, byrow = TRUE, nrow = n_sample, ncol = 4)
  
  dat_X = cbind(
    dat_X_X1V1X2V2[, c(1, 3)], dat_X_X3,
    dat_X_X1V1X2V2[, c(2, 4)], dat_X_V3
  )
  colnames(dat_X) = c(paste0('X', 1:3), paste0('V', 1:3))
  
  return(dat_X)
  
}

para_PS_fn = function(X, level = 'str') {
  # input 
  #   X: covariate vector X (X1-X3, V1-V3)
  #   level: string, 'str' (strong) or 'mod' (moderate)
  #     for the association between X and W
  
  if (level == 'str') {
    para_PS_beta_0 = 0
    para_PS_beta_X = c(0.6, -0.6, 0.6)
  } else if (level == 'mod') {
    para_PS_beta_0 = 0
    para_PS_beta_X = c(0.3, -0.3, 0.3)
  }
  
  dat_PS = plogis(para_PS_beta_0 + sum(X[1:3] * para_PS_beta_X))
  
  return(dat_PS)
  
}

para_Y_fn = function(W, X, level = 'str') {
  # input: 
  #   W: binary, treatment indicator
  #   X: covariate vector X (X1-X3, V1-V3)
  #   level: string, 'str' (strong) or 'mod' (moderate) or 'no'
  #     for the association between V and Y
  
  para_Y_TE = 2
  para_Y_beta_0 = 0
  para_Y_beta_X = c(-1, 1, -1)
  if (level == 'str') {
    para_Y_beta_V = c(-1, 1, 1)
  } else if (level == 'mod') {
    para_Y_beta_V = c(-0.5, 0.5, 0.5)
  } else if (level == 'no') {
    para_Y_beta_V = c(0, 0, 0)
  }
  
  Y = para_Y_beta_0 + para_Y_TE * W + 
    sum(para_Y_beta_X * X[1:3]) + sum(para_Y_beta_V * X[4:6])
  
  return(Y)
  
}

para_Y_sigma_error = 1

# Simulate one dataset ----------------------------------------------------

set.seed(20200512) 

para_level_X_PS = 'str'
para_level_X_Y = 'str'

dat_X = para_X_fn(n_sample)

dat_PS = apply(dat_X, 1, para_PS_fn, level = para_level_X_PS)
dat_W = rbinom(n_sample, size = 1, prob = dat_PS)

dat_Y0 = apply(dat_X, 1, para_Y_fn, W = 0, level = para_level_X_Y)
dat_Y1 = apply(dat_X, 1, para_Y_fn, W = 1, level = para_level_X_Y)
dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(n_sample, sd = para_Y_sigma_error)
dat_Y = cbind(dat_Y0, dat_Y1, dat_Y_obs)
colnames(dat_Y) = c('Y0', 'Y1', 'Y')

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)