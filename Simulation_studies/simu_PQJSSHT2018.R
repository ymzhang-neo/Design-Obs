# PQJSSHT2018
# Some methods for heterogeneous treatment effect estimation in high dimensions
# Statistics in Medicine

set.seed(20201103)


# 0 Preparation -----------------------------------------------------------

# Settings of the scenarios 

scenario_n_sample_all = c(200, 200, 300, 300, 400, 400, 1000, 1000,
                          200, 200, 300, 300, 400, 400, 1000, 1000)

scenario_n_cov_all = c(400, 400, 300, 300, 200, 200, 100, 100,
                       400, 400, 300, 300, 200, 200, 100, 100)

scenario_mean_fn_id = c(8, 5, 4, 7, 3, 1, 2, 6,
                        8, 5, 4, 7, 3, 1, 2, 6)

scenario_ite_fn_id = c(1:8, 1:8)

scenario_Y_var = c(1, 1/4, 1, 1/4, 1, 1, 4, 4, 
                   1, 1/4, 1, 1/4, 1, 1, 4, 4)

# Building-block functions

fn_f1 = function(x) {0}

fn_f2 = function(x) {5*(x[1] > 1) - 5}

fn_f3 = function(x) {2*x[1] - 4}

fn_f4 = function(x) {
  x[2]*x[4]*x[6] + 2*x[2]*x[4]*(1-x[6]) + 
    3*x[2]*(1-x[4])*x[6] + 4*x[2]*(1-x[4])*(1-x[6]) + 
    5*(1-x[2])*x[4]*x[6] + 6*(1-x[2])*x[4]*(1-x[6]) + 
    7*(1-x[2])*(1-x[4])*x[6] + 8*(1-x[2])*(1-x[4])*(1-x[6])
}

fn_f5 = function(x) {x[1]+x[3]+x[5]+x[7]+x[8]+x[9]-2}

fn_f6 = function(x) {
  4*(x[1] > 1)*(x[3] > 0) + 4*(x[5] > 1)*(x[7] > 0) + 2*x[8]*x[9]
}

fn_f7 = function(x) {
  (x[1]^2 + x[2] + x[3]^2 + x[4] + x[5]^2 + x[6] + x[7]^2 + x[8] + x[9]^2 
   - 11) / 2
}

fn_f8 = function(x) {(fn_f4(x) + fn_f5(x)) / sqrt(2)}

fn_all = list(fn_f1, fn_f2, fn_f3, fn_f4, fn_f5, fn_f6, fn_f7, fn_f8)


# 1 Generate data ---------------------------------------------------------

# Choose the scenario and the parameters

i_scenario = 1 # 1:16

n_sample = scenario_n_sample_all[i_scenario]

n_cov = scenario_n_cov_all[i_scenario]

mean_fn_id = scenario_mean_fn_id[i_scenario]

ite_fn_id = scenario_ite_fn_id[i_scenario]

Y_var = scenario_Y_var[i_scenario]

# Generate covariates

dat_X = as.data.frame(matrix(NA, nrow = n_sample, ncol = n_cov))

colnames(dat_X) = paste0('X', 1:n_cov)

dat_X[, seq(from=1, to=n_cov, by=2)] = matrix(
  rnorm(n_cov * n_sample / 2), nrow = n_sample
)

dat_X[, seq(from=2, to=n_cov, by=2)] = matrix(
  rbinom(n_cov * n_sample / 2, size = 1, prob = 0.5), nrow = n_sample
)

# Generate ITE 

dat_Y = as.data.frame(matrix(NA, nrow = n_sample, ncol = 7))

colnames(dat_Y) = c('mu_mean', 'ITE', 'mu0', 'mu1', 'Y0', 'Y1', 'Y_obs')

dat_Y$mu_mean = apply(dat_X, 1, fn_all[[mean_fn_id]])

dat_Y$ITE = apply(dat_X, 1, fn_all[[ite_fn_id]])

dat_Y$mu0 = dat_Y$mu_mean - (dat_Y$ITE / 2)

dat_Y$mu1 = dat_Y$mu_mean + (dat_Y$ITE / 2)

dat_Y$Y0 = dat_Y$mu0 + rnorm(n_sample, sd = sqrt(Y_var))

dat_Y$Y1 = dat_Y$mu1 + rnorm(n_sample, sd = sqrt(Y_var))

# Generate PS and W

dat_PS = as.data.frame(matrix(NA, nrow = n_sample, ncol = 2))

colnames(dat_PS) = c('PS_logit', 'PS')

if (i_scenario <= 8) {
  dat_PS$PS_logit = rep(0, n_sample)
} else {
  dat_PS$PS_logit = dat_Y$mean - (dat_Y$ITE / 2)
}

dat_PS$PS = plogis(dat_PS$PS_logit)

dat_W = rbinom(n_sample, size = 1, prob = dat_PS_prob)

# Generate Y_obs

dat_Y$Y_obs = dat_W * dat_Y$Y1 + (1 - dat_W) * dat_Y$Y0

# Combine all together 

dat = data.frame(dat_X, dat_PS, W = dat_W, dat_Y)

