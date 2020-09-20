# Rudolph, Colson, Stuart, Ahern 2016
# Optimally combining propensity score subclasses
# Statistics in Medicine

set.seed(20200919)

# Covariates --------------------------------------------------------------

n_stra = 5 # 10, 30

n_sample = ifelse(n_stra == 30, 5000, 2000)

dat_X = as.data.frame(matrix(NA, nrow = n_sample, ncol = 4))
colnames(dat_X) = paste0('X', 1:4)
dat_X$X1 = runif(n_sample, min = 0.02, max = 0.7)
dat_X$X2 = rnorm(n_sample, mean = 0.2 + 0.125 * dat_X$X1, sd = 1)
dat_X$X3 = rnorm(n_sample, mean = -2, sd = 0.7)
dat_X$X4 = rbinom(n_sample, size = 1, prob = 0.4)

dat_survey = as.data.frame(matrix(NA, nrow = n_sample, ncol = 3))
colnames(dat_survey) = c('incl_logit', 'incl_prob', 'incl_ind')
dat_survey$incl_logit = dat_X$X1 - dat_X$X3 + dat_X$X4
dat_survey$incl_prob = plogis(dat_survey$incl_logit)
dat_survey$incl_ind = rbinom(n_sample, size = 1, prob = dat_survey$incl_prob)

dat_Y = as.data.frame(matrix(NA, nrow = n_sample, ncol = 5))
colnames(dat_Y) = c('mu0', 'mu1', 'Y0', 'Y1', 'Y_obs')

# DGM 1 -------------------------------------------------------------------
# Base scenario: constant treatment effect; positivity assumption

dat_PS_logit = -0.5 + dat_X$X1 + 0.1 * dat_X$X1^2 - 0.5 * dat_X$X2 + 
  0.5 * dat_X$X1 * dat_X$X2
dat_PS_prob = plogis(dat_PS_logit)
dat_W = rbinom(n_sample, size = 1, dat_PS_prob)

dat_Y$mu0 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2
dat_Y$mu1 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2 + 2


# DGM 2 -------------------------------------------------------------------
# Constant treatment effect; violation of positivity assumption

dat_PS_logit = -0.3 + dat_X$X1 - 1.5 * dat_X$X1^2 - 1.5 * dat_X$X2 + 
  1.5 * dat_X$X1 * dat_X$X2
dat_PS_prob = plogis(dat_PS_logit)
dat_W = rbinom(n_sample, size = 1, dat_PS_prob)

dat_Y$mu0 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2
dat_Y$mu1 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2 + 2


# DGM 3 -------------------------------------------------------------------
# Heterogeneous treatment effect; positivity assumption

dat_PS_logit = -0.5 + dat_X$X1 + 0.1 * dat_X$X1^2 - 0.5 * dat_X$X2 + 
  0.5 * dat_X$X1 * dat_X$X2
dat_PS_prob = plogis(dat_PS_logit)
dat_W = rbinom(n_sample, size = 1, dat_PS_prob)

dat_Y$mu0 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2
dat_Y$mu1 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2 + 
  1.5 + 2 * dat_X$X1 + 2 * dat_X$X1^2 - dat_X$X2


# DGM 4 -------------------------------------------------------------------
# Heterogeneous treatment effect; violation of positivity assumption
# [N] This is an additional scenario created based on DGM 2 and 3.

dat_PS_logit = -0.3 + dat_X$X1 - 1.5 * dat_X$X1^2 - 1.5 * dat_X$X2 + 
  1.5 * dat_X$X1 * dat_X$X2
dat_PS_prob = plogis(dat_PS_logit)
dat_W = rbinom(n_sample, size = 1, dat_PS_prob)

dat_Y$mu0 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2
dat_Y$mu1 = -0.5 + 3 * dat_X$X1 + 3 * dat_X$X1^2 - 2 * dat_X$X2 + 
  1.5 + 2 * dat_X$X1 + 2 * dat_X$X1^2 - dat_X$X2


# Combining ---------------------------------------------------------------

dat_Y$Y0 = dat_Y$mu0 + rnorm(n_sample, mean = 0, sd = 1)
dat_Y$Y1 = dat_Y$mu1 + rnorm(n_sample, mean = 0, sd = 1)
dat_Y$Y_obs = dat_W * dat_Y$Y1 + (1 - dat_W) * dat_Y$Y0

dat = data.frame(dat_X, PS_logit = dat_PS_logit, PS = dat_PS_prob, W = dat_W, 
                 dat_Y)


# Understanding -----------------------------------------------------------

# PS and covariate balance

plot_X = 'X2'
plot_xlim = range(dat[, plot_X])
par(mfrow = c(2, 1))
hist(dat[dat$W == 0, plot_X], xlim = plot_xlim, xlab = plot_X, main = 'W = 0')
hist(dat[dat$W == 1, plot_X], xlim = plot_xlim, xlab = plot_X, main = 'W = 1')
par(mfrow = c(1, 1))

# Survey membership
# The authors considered survey data where some subjects are not observed. 
# However, the distributions do not look much different, as only 10% subjects
#   are discarded.

plot_X = 'X2'
plot_xlim = range(dat[, plot_X])
par(mfrow = c(2, 2))
hist(dat[dat$W == 0, plot_X], xlim = plot_xlim, xlab = plot_X, main = 'W = 0')
hist(dat[(dat$W == 0) & (dat_survey$incl_ind == 1), plot_X], 
     xlim = plot_xlim, xlab = plot_X, main = 'W = 0, Observed')
hist(dat[dat$W == 1, plot_X], xlim = plot_xlim, xlab = plot_X, main = 'W = 1')
hist(dat[(dat$W == 0) & (dat_survey$incl_ind == 1), plot_X], 
     xlim = plot_xlim, xlab = plot_X, main = 'W = 1, Observed')
par(mfrow = c(1, 1))

table(dat$W, dat_survey$incl_ind)
