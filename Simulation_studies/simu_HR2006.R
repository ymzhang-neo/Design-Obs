# Hill, Reiter (2006)
# Interval estimation of treatment effects using propensity score matching
# Statistics in Medicine

# Parameters --------------------------------------------------------------

n_simu = 1000
n_boot = 250

para_Y_mean_fn = list(
  linear_additive = function(X, W) {
    4 * W + X[1] + 2 * X[2]
  },
  nonlinear_additive = function(X, W) {
    4 * W + X[1] * X[2] + 0.5 * X[1]^2 + 1.5 * X[2]*2
  }, 
  nonlinear_nonadditive = function(X, W) {
    if (W == 0) {
      X[1] + 2 * X[2]
    } else if (W == 1) {
      4 + X[1] * X[2] + 0.5 * X[1]^2 + 1.5 * X[2]*2
    }
  }
)

para_Y_sigma_error = 1

#### Scenario 1: strong overlap ####

n_sample_W1 = 150
n_sample_W0 = 350

dat_X = matrix(rnorm((n_sample_W1 + n_sample_W0) * 2, mean = 1, sd = 1), ncol = 2)
colnames(dat_X) = c('X1', 'X2')

dat_W = c(rep(1, n_sample_W1), rep(0, n_sample_W0))
dat_PS = rep(NA, length(dat_W))

i_Y_mean = 1 # 2, 3

dat_Y0 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 0)
dat_Y1 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 1)
dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(n_sample_W0 + n_sample_W1, 
        mean = 0, sd = para_Y_sigma_error)
# dat_Y_obs = mapply(para_Y_mean_fn[[1]], as.data.frame(t(dat_X)), dat_W)
dat_Y = cbind(Y0 = dat_Y0, Y1 = dat_Y1, Y = dat_Y_obs)

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)

#### Scenario 2: moderate overlap ####

n_sample_W1 = 150
n_sample_W0 = 150
n_sample_W0_distractor = 200

dat_X = rbind(
  matrix(rnorm((n_sample_W1 + n_sample_W0) * 2, mean = 1, sd = 1), ncol = 2), 
  matrix(rnorm(n_sample_W0_distractor * 2, mean = 3, sd = 1), ncol = 2)
)
colnames(dat_X) = c('X1', 'X2')

dat_W = c(rep(1, n_sample_W1), rep(0, n_sample_W0 + n_sample_W0_distractor))
dat_PS = rep(NA, length(dat_W))

i_Y_mean = 1 # 2, 3

dat_Y0 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 0)
dat_Y1 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 1)
dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(n_sample_W0 + n_sample_W1 + n_sample_W0_distractor, 
        mean = 0, sd = para_Y_sigma_error)
dat_Y = cbind(dat_Y0, dat_Y1, dat_Y_obs)
colnames(dat_Y) = c('Y0', 'Y1', 'Y')

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)

#### Scenario 3: weak overlap ####

n_sample_W1 = 150
n_sample_W0 = 50
n_sample_W0_distractor = 300

dat_X = rbind(
  matrix(rnorm((n_sample_W1 + n_sample_W0) * 2, mean = 1, sd = 1), ncol = 2), 
  matrix(rnorm(n_sample_W0_distractor * 2, mean = 3, sd = 1), ncol = 2)
)
colnames(dat_X) = c('X1', 'X2')

dat_W = c(rep(1, n_sample_W1), rep(0, n_sample_W0 + n_sample_W0_distractor))
dat_PS = rep(NA, length(dat_W))

i_Y_mean = 1 # 2, 3

dat_Y0 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 0)
dat_Y1 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 1)
dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(n_sample_W0 + n_sample_W1 + n_sample_W0_distractor, 
        mean = 0, sd = para_Y_sigma_error)
dat_Y = cbind(Y0 = dat_Y0, Y1 = dat_Y1, Y = dat_Y_obs)  

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)

#### Scenario 4: uneven overlap ####

n_sample = 500 

dat_X = matrix(runif(n_sample * 2, min = -1, max = 3), ncol = 2)
colnames(dat_X) = c('X1', 'X2')

dat_PS_category = cut(rowSums(dat_X), breaks = c(-2, 0, 1, 3, 4, 6), labels = FALSE)
para_PS = c(0.05, 0.10, 0.15, 0.65, 0.75)
dat_PS = para_PS[dat_PS_category]
dat_W = rbinom(n_sample, size = 1, prob = dat_PS)

i_Y_mean = 2 # 2, 3

dat_Y0 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 0)
dat_Y1 = apply(dat_X, 1, para_Y_mean_fn[[i_Y_mean]], W = 1)
dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(n_sample, 
        mean = 0, sd = para_Y_sigma_error)
dat_Y = cbind(Y0 = dat_Y0, Y1 = dat_Y1, Y = dat_Y_obs)

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)

# Estimate without design -------------------------------------------------

mean(dat[dat$W == 1, 'Y']) - mean(dat[dat$W == 0, 'Y'])

# Notes -------------------------------------------------------------------
# - About potential outcomes 
#     Y_obs = (1-W)*Y0 + W*Y1
#     Y_mis = W*Y0 + (1-W)*Y1
#     Y0 = (1-W)*Y_obs * W*Y_mis
#     Y1 = W*Y_obs + (1-W)*Y_mis
#   (No random error assumed.)