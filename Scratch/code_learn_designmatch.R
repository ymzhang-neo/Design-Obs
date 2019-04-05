# code_learn_designmatch

# https://cran.rstudio.com/web/packages/designmatch/designmatch.pdf

library(designmatch)


# Design: bmatch() --------------------------------------------------------

data(lalonde)
lalonde = lalonde[order(lalonde$treatment, decreasing = TRUE), ]
attach(lalonde)

#### bmatch() - cardinality matching ####

t_ind = treatment

dist_mat = NULL # distance matrix
subset_weight = 1 # subset matching weight

# moment balance
mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
mom = list(covs = mom_covs, tols = mom_tols)

# fine balance
fine_covs = cbind(black, hispanic, married, nodegree)
fine = list(covs = fine_covs)

# exact matching 
exact_covs = cbind(black)
exact = list(covs = exact_covs)

# solver options
t_max = 60*5
solver = 'glpk'
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate, 
              round_cplex = 0,  trace = 0)

# match!
out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight, 
             mom = mom, fine = fine, exact = exact, solver = solver)
str(out)

# assess mean balance
meantab(mom_covs, t_ind, out$t_id, out$c_id)

# assess fine balance
for (i in 1:ncol(fine_covs)) {
  print(finetab(fine_covs[, i], out$t_id, out$c_id))
}

# assess exact matching balance
table(exact_covs[out$t_id] == exact_covs[out$c_id])


# Covariate balance -------------------------------------------------------

# absstddif()
# [Q] meaning of std_dif=.
data(lalonde)
attach(lalonde)
t_ind = treatment
mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
mom_tols = absstddif(X_mat = mom_covs, t_ind = t_ind, std_dif = 1)
