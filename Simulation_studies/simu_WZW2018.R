# Simulation: WZW2018
# Wu, Zeng, Wang (2018)
# Matched learning for optimizing individualized treatment strategies using electronic health records
# Journal of the American Statistical Association

# Overview ----------------------------------------------------------------
# Notes -------------------------------------------------------------------
# - In Designs 3 and 4,  it is not clear how the unmeasured covariate X 
#   should be simulated.


para_N = 100 # 200, 500, 1000

# Design 1 ----------------------------------------------------------------

fn_s1 = function(x, W) {
  W * (x[1] - x[2]) + 6 * sign(x[1]) + 2 * x[3] - x[4] + 
  rnorm(1)
}

fn_ps1 = function(x, type = 'logit') {
  ps_logit = 1 + 2 * x[1] + x[2]
  ps_prob = exp(ps_logit) / (1 + exp(ps_logit))
  if (type == 'logit') {
    ps = ps_logit
  } else if (type == 'prob') {
    ps = ps_prob
  }
  return(ps)
}

dat_X = matrix(rnorm(para_N * 4), ncol = 4)

dat_W_prob = apply(dat_X, 1, fn_ps1, type = 'prob')
dat_W = rbinom(para_N, size = 1, prob = dat_W_prob)

dat_Y = sapply(1:para_N, function(oneID) {fn_s1(x = dat_X[oneID, ], W = dat_W[oneID])})


# Design 2 ----------------------------------------------------------------

fn_s2 = function(x, W) {
  W * (x[1]*2 + x[2] - 1) + 1 + 2 * x[1] + x[2] + 0.5 * x[3] + 6 * sign(x[1]) + 
  rnorm(1)
}

fn_ps2 = function(x, type = 'logit') {
  ps_logit = 1 + exp(x[2])
  ps_prob = exp(ps_logit) / (1 + exp(ps_logit))
  if (type == 'logit') {
    ps = ps_logit
  } else if (type == 'prob') {
    ps = ps_prob
  }
  return(ps)
}

dat_X = matrix(rnorm(para_N * 3), ncol = 3)

dat_W_prob = apply(dat_X, 1, fn_ps2, type = 'prob')
dat_W = rbinom(para_N, size = 1, prob = dat_W_prob)

dat_Y = sapply(1:para_N, function(oneID) {fn_s2(x = dat_X[oneID, ], W = dat_W[oneID])})


# Design 3 ----------------------------------------------------------------

fn_s3_wNe = function(x, x_unmeasured) {
  -(x[1] - x[2] + x_unmeasured) + 6 * sign(x[1]) + 2 * x[3] - x[4]
}
fn_s3_wPo = function(x, x_unmeasured) {
  (x[1] - x[2] + x_unmeasured) + 6 * sign(x[1]) + 2 * x[3] - x[4]
}

fn_ps3 = function(x, x_unmeasured, type = 'logit') {
  ps_logit = 1 + x[1] + 2 * x_unmeasured + 
    fn_s3_wNe(x = x, x_unmeasured = x_unmeasured) - 
    fn_s3_wPo(x = x, x_unmeasured = x_unmeasured) 
  ps_prob = exp(ps_logit) / (1 + exp(ps_logit))
  if (type == 'logit') {
    ps = ps_logit
  } else if (type == 'prob') {
    ps = ps_prob
  }
  return(ps)
}

dat_X = matrix(rnorm(para_N * 4), ncol = 4)

dat_X_unmeasured = rnorm(para_N)

dat_W_prob = sapply(1:para_N, function(oneID) {
  fn_ps3(x = dat_X[oneID, ], x_unmeasured = dat_X_unmeasured[oneID], type = 'prob')
})
dat_W = rbinom(para_N, size = 1, prob = dat_W_prob)

dat_Y = sapply(1:para_N, function(oneID) {
  if (dat_W[oneID] == 0) {
    oneY = fn_s3_wNe(x = dat_X[oneID, ], x_unmeasured = dat_X_unmeasured[oneID])
  } else if (dat_W[oneID] == 1) {
    oneY = fn_s3_wPo(x = dat_X[oneID, ], x_unmeasured = dat_X_unmeasured[oneID])
  }
}) + rnorm(para_N)


# Design 4 ----------------------------------------------------------------

fn_s4_wNe = function(x, x_unmeasured) {
  -(x[1]^2 + x[2] + x_unmeasured - 1) + 6 * sign(x[1]) + 1 + 2 * x[1] + x[2]
}
fn_s4_wPo = function(x, x_unmeasured) {
  (x[1]^2 + x[2] + x_unmeasured - 1) + 6 * sign(x[1]) + 1 + 2 * x[1] + x[2]
}

fn_ps4 = function(x, x_unmeasured, type = 'logit') {
  ps_logit = 1 + x[1] + 2 * x_unmeasured + 
    fn_s4_wNe(x = x, x_unmeasured = x_unmeasured) - 
    fn_s4_wPo(x = x, x_unmeasured = x_unmeasured) 
  ps_prob = exp(ps_logit) / (1 + exp(ps_logit))
  if (type == 'logit') {
    ps = ps_logit
  } else if (type == 'prob') {
    ps = ps_prob
  }
  return(ps)
}

dat_X = matrix(rnorm(para_N * 4), ncol = 4)

dat_X_unmeasured = rnorm(para_N)

dat_W_prob = sapply(1:para_N, function(oneID) {
  fn_ps4(x = dat_X[oneID, ], x_unmeasured = dat_X_unmeasured[oneID], type = 'prob')
})
dat_W = rbinom(para_N, size = 1, prob = dat_W_prob)

dat_Y = sapply(1:para_N, function(oneID) {
  if (dat_W[oneID] == 0) {
    oneY = fn_s4_wNe(x = dat_X[oneID, ], x_unmeasured = dat_X_unmeasured[oneID])
  } else if (dat_W[oneID] == 1) {
    oneY = fn_s4_wPo(x = dat_X[oneID, ], x_unmeasured = dat_X_unmeasured[oneID])
  }
}) + rnorm(para_N)
