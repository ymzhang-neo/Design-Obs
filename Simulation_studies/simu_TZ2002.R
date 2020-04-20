# Tu, Zhou (2003)
# A bootstrap confidence interval procedure for the treatment effect using propensity score subclassification
# Health Services & Outcomes Research Methodology

# Parameters --------------------------------------------------------------

size_W0 = c(500, 1000, 2000)
size_W1 = size_W0

para_df = data.frame(
  # PS model
  d = 0.5,
  sigma_X1 = rep(c(1, 1.25), each = 4), 
  p_X2_W0 = rep(c(0.3, 0.4), each = 4),
  p_X3_W0 = rep(c(0.7, 0.8), each = 4),
  p_X2_W1 = rep(c(0.8, 0.7), each = 4),
  p_X3_W1 = rep(c(0.4, 0.3), each = 4),
  # Y model
  TE = rep(c(0, 0.5, 1, 2), 2), 
  beta1 = 0.5, 
  beta2 = 0.4, 
  beta3 = 0.4,
  sigma_error = 1.0
)

n_simu = 1000
n_boot = 2000

  
# Simulate one dataset ----------------------------------------------------

set.seed(20200407+1832)

i_para = 1 # from 1:8
i_size = 1 # from 1:3

# PS and W

dat_PS = rep(NA, sum(size_W0[i_size], size_W1[i_size]))
dat_W = rep(c(0, 1), c(size_W0[i_size], size_W1[i_size]))

# X

dat_X1 = c(
  rnorm(size_W0[i_size], mean = 0, 
        sd = para_df[i_para, 'sigma_X1']), 
  rnorm(size_W1[i_size], mean = para_df[i_para, 'd'], 
        sd = para_df[i_para, 'sigma_X1'])
)
dat_X2 = c(
  rbinom(size_W0[i_size], size = 1, prob = para_df[i_para, 'p_X2_W0']), 
  rbinom(size_W1[i_size], size = 1, prob = para_df[i_para, 'p_X2_W1'])
)
dat_X3 = c(
  rbinom(size_W0[i_size], size = 1, prob = para_df[i_para, 'p_X3_W0']), 
  rbinom(size_W1[i_size], size = 1, prob = para_df[i_para, 'p_X3_W1'])
)
dat_X = cbind(dat_X1, dat_X2, dat_X3)
colnames(dat_X) = c('X1', 'X2', 'X3')

# Y 

dat_Y0 = as.vector(dat_X %*% t(para_df[i_para, c('beta1', 'beta2', 'beta3')])) 
dat_Y1 = para_df[i_para, 'TE'] + dat_Y0
dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(size_W0[i_size] + size_W1[i_size], sd = para_df[i_para, 'sigma_error'])
dat_Y = cbind(dat_Y0, dat_Y1, dat_Y_obs)
colnames(dat_Y) = c('Y0', 'Y1', 'Y')

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)

# Estimate without design

D0_estN_est = mean(dat[dat$W == 1, 'Y']) - mean(dat[dat$W == 0, 'Y'])

# Notes -------------------------------------------------------------------
# - The authors did not generate true PS.
# - About potential outcomes (not used here)
#     Y_obs = (1-W)*Y0 + W*Y1
#     Y_mis = W*Y0 + (1-W)*Y1
#     Y0 = (1-W)*Y_obs * W*Y_mis
#     Y1 = W*Y_obs + (1-W)*Y_mis
#   (No random error assumed.)