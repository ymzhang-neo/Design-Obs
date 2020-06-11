# Simulation: WA2018
# Wager, Athey (2018)
# Estimation and inference of heterogeneous treatment effects using random forests
# Journal of the American Statistical Association

# Overview ----------------------------------------------------------------

## The authors proposed four designs. In each design, they specified 
##   - d: covariate dimension
##   - N: total sample size
##   - e: propensity score
##   - m: mean-effect function (input: X)
##   - tau: ITE function (input: X)
## With the setup, the mean functions for the potential outcomes can be 
##   recovered: 
##   - E[Y(0) | X] = fn_m(X) - fn_tau(X) / 2
##   - E[Y(1) | X] = fn_m(X) + fn_tau(X) / 2
##   The specific potential outcome values have noise from N(0, 1) added.


# Parameters --------------------------------------------------------------

#### Design 1 ####

# para_d1_all = c(2, 5, 10, 15, 20, 30)
# para_N1 = 500 

fn_e1 = function(x) {(1 + dbeta(x[1], shape1 = 2, shape2 = 4)) / 4}

fn_m1 = function(x) {2 * x[1] - 1}

fn_tau1 = function(x) {0}

#### Design 2 ####

# para_d2_all = c(2, 3, 4, 5, 6, 8)
# para_N2 = 5000

fn_e2 = function(x) {0.5}

fn_m2 = function(x) {0}

fn_tau2 = function(x) {
  (1 + 1 / (1 + exp(-20 * (x[1] - 1 / 3)))) * 
  (1 + 1 / (1 + exp(-20 * (x[2] - 1 / 3))))
} 

#### Design 3 ####

# para_d3_all = c(2, 3, 4, 5, 6, 8)
# para_N3 = 10000

fn_e3 = function(x) {0.5}

fn_m3 = function(x) {0}

fn_tau3 = function(x) {
  (1 + 2 / (1 + exp(-12 * (x[1] - 1 / 2)))) * 
  (1 + 2 / (1 + exp(-12 * (x[2] - 1 / 2))))
} 

#### Design 4 #### (in appendix, more work needed)

# para_d4_all = c(6, 12)
# para_N4 = 5000

fn_e4 = function(x) {0.5}
fn_m4 = function(x) {0}

# Treatment effect, x = input vector of length d
# q = number of covariates with signal, 2, 4, 6
fn_varsigma4 = function(x) {(1 / (1 + exp(-12 * (x - 0.5))) - 0.5)}
fn_tau4 = function(x, q) {sum(sapply(x[1:q], fn_varsigma4)) * 4 / q}

#### Misc. ####

fn_design_all = list(
  list(e = fn_e1, m = fn_m1, tau = fn_tau1),
  list(e = fn_e2, m = fn_m2, tau = fn_tau2),
  list(e = fn_e3, m = fn_m3, tau = fn_tau3),
  list(e = fn_e4, m = fn_m4, tau = fn_tau4, fn_varsigma4)
)

para_Y_sigma_error = 1

# Dimension and Size 
para_d_all = list(
  c(2, 5, 10, 15, 20, 30), 
  c(2, 3, 4, 5, 6, 8),
  c(2, 3, 4, 5, 6, 8),
  c(6, 12)
)
para_n_sample_all = c(500, 5000, 10000, 5000)


# Get data ----------------------------------------------------------------

set.seed(1)

i_design = 1

d = para_d_all[[i_design]][1]

n_sample = para_n_sample_all[i_design]

dat_X = matrix(runif(n_sample * d), ncol = d)
colnames(dat_X) = paste0('X', 1:d)

dat_PS = apply(dat_X, 1, fn_design_all[[i_design]][['e']])

dat_W = rbinom(n_sample, size = 1, prob = dat_PS)

dat_Y = data.frame(
  m = apply(dat_X, 1, fn_design_all[[i_design]][['m']]), 
  tau = apply(dat_X, 1, fn_design_all[[i_design]][['tau']])
)
dat_Y$Y0_mean = dat_Y$m - dat_Y$tau / 2
dat_Y$Y1_mean = dat_Y$m + dat_Y$tau / 2

dat_Y$Y_obs_mean = dat_W * dat_Y$Y1_mean + (1 - dat_W) * dat_Y$Y0_mean
dat_Y$Y_mis_mean = (1 - dat_W) * dat_Y$Y1_mean + dat_W * dat_Y$Y0_mean
dat_Y$Y_obs = dat_Y$Y_obs_mean + rnorm(n_sample, sd = para_Y_sigma_error)
dat_Y$Y_mis = dat_Y$Y_mis_mean + rnorm(n_sample, sd = para_Y_sigma_error)
dat_Y$Y0 = (1 - dat_W) * dat_Y$Y_obs + dat_W * dat_Y$Y_mis
dat_Y$Y1 = dat_W * dat_Y$Y_obs + (1 - dat_W) * dat_Y$Y_mis

dat_true = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)

dat = data.frame(dat_X, W = dat_W, Y = dat_Y$Y_obs)

## Notes
## - The current code does not allow easy specification of q in Design 4.


# Visualization -----------------------------------------------------------

par(mfrow = c(2, 1))

xLim = range(dat_Y)
hist(dat_Y[dat_W == 0], xlab = 'Y', xlim = xLim, main = 'W == 0')
hist(dat_Y[dat_W == 1], xlab = 'Y', xlim = xLim, main = 'W == 1')

par(mfrow = c(1, 1))