# learn_optmatch
# Code for learning to use the package optmatch & RItools (HFP2013)

library(optmatch)
library(RItools)

data('nuclearplants')
head(nuclearplants)

# pt is not the treatment considered in this example, yet.
table(nuclearplants$pt)
# with(nuclearplants, table(pt))

nuke.nopt = subset(nuclearplants, pt == 0)
table(nuke.nopt$pr)

pm = pairmatch(match_on(pr ~ cap, data = nuke.nopt))
summary(pm)
matched(pm)

cap.dist = match_on(pr ~ cap, data = nuke.nopt)

cap.dist[, 1:3]




# Scratch -----------------------------------------------------------------

dat_forPairm = cbind(W = dat_forExact_ID$W, dat_forExact)
temp.pairm = pairmatch(match_on(W ~ ., data = dat_forPairm))
summary(temp.pairm)

temp.mydat = data.frame(
  W = c(0, 0, 1, 1, 1),
  PS_link = c(0, 1, 2, 3, -1)
)

outer(temp.mydat[temp.mydat$W == 0, 'PS_link'], 
      temp.mydat[temp.mydat$W == 1, 'PS_link'], 
      FUN = function(x, y){abs(x - y)})
