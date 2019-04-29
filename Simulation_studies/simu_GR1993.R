# Simulation: GR1993
# Gu, Rosenbaum (1993)
# Comparison of Multivariate Matching Methods: Structures, Distances, and Algorithms
# Journal of Computational and Graphical Statistics


# Overview ----------------------------------------------------------------

# The article designed simulation settings to investigate three aspects of 
# matching procedures: algorithm, structure, and distance.

set.seed(1)

para_r = 1 # 2, 3, 6 # Ratio of #control:#treated
para_N_W1 = 50
para_N_W0 = para_N_W1 * para_r
para_N = para_N_W0 + para_N_W1

para_p = 2 # 5, 20

para_b = 0.5 # 1.0 # bias
para_mu_W1 = c(para_b, rep(0, para_p-1))
para_mu_W0 = rep(0, para_p)

dat_X = rbind(
  t(para_mu_W1 + matrix(rnorm(para_N_W1*para_p), nrow = para_p)),
  t(para_mu_W0 + matrix(rnorm(para_N_W0*para_p), nrow = para_p))
)
colnames(dat_X) = paste0('X', 1:para_p)

dat_W = c(rep(1, para_N_W1), rep(0, para_N_W0))

