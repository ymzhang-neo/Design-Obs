# Austin, Small (2014)
# The use of bootstrapping when using propensity-score matching without replacemnet: a simulation study
# Statistics in Medicine 

# Parameters --------------------------------------------------------------

n_sample = 5000
n_simu = 1000
n_boot = 1000

para_X_p = 10

para_beta = c(log(1.25), log(1.5), log(1.75), log(2))
names(para_beta) = c('W', 'M', 'S', 'VS') # weak, moderate, strong, very strong

para_PS_beta_X = c(rep(para_beta[-4], 2), para_beta['VS'])
para_PS_W1_prop = c(0.05, 0.1, 0.2, 0.25)
para_PS_beta0 = log(para_PS_W1_prop / (1 - para_PS_W1_prop))

para_Y_TE = 1
para_Y_beta_X = c(para_beta, para_beta[-4])
para_Y_sigma_error = 3

# Simulate one dataset ----------------------------------------------------

set.seed(20200409 + 2154)

dat_X = matrix(rnorm(n_sample * para_X_p), ncol = para_X_p)
colnames(dat_X) = sprintf("X%02d", 1:para_X_p)

i_W1_prop = 1 # 1, 2, 3, 4

dat_PS_logit = para_PS_beta0[i_W1_prop] + 
  dat_X[, 1:7] %*% para_PS_beta_X
dat_PS = plogis(dat_PS_logit)
dat_W = rbinom(n_sample, size = 1, prob = dat_PS)

dat_Y0 = dat_X[, 4:10] %*% para_Y_beta_X
dat_Y1 = dat_Y0 + para_Y_TE
dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(n_sample, sd = para_Y_sigma_error)
dat_Y = cbind(dat_Y0, dat_Y1, dat_Y_obs)
colnames(dat_Y) = c('Y0', 'Y1', 'Y')

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)

# Notes -------------------------------------------------------------------
# - Outcomes
#   The authors also considered binary and survival outcomes. 
# - Matching algorithms
#   - greedy nearest-neighbor matching 
#   - caliper matching
#   - optimal matching


