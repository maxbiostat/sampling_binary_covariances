source("aux.r")
library(hitandrun)
library(bindata)
###########################

P <- 8
mprobs <- runif(P)

make_constraints_har_pairwiseOnly <- function(pairwise.constraints,
                                              ncores = 6) {
  ### Pairwise stuff
  
  prw.tab <- pairwise.constraints$pairwise_ctr
  ind.mat <- pairwise.constraints$ind_mat
  M <- nrow(prw.tab)
  
  l.prw <- parallel::mclapply(1:M, function(t) {
    pp <- prw.tab[t,]$par_index
    lowerBoundConstraint(n = M, i = pp, x = prw.tab[t,]$lb)
  }, mc.cores = ncores)
  
  u.prw <- parallel::mclapply(1:M, function(t) {
    pp <- prw.tab[t,]$par_index
    upperBoundConstraint(n = M, i = pp, x = prw.tab[t,]$ub)
  }, mc.cores = ncores)
  
  
  
  the.ctrs <- Reduce(function(x, y)
    mergeConstraints(x, y),
    c(l.prw, u.prw))
  return(the.ctrs)
}

pairwise <- get_pairwise_constraints(ps = mprobs)

all.constraints <- make_constraints_har_pairwiseOnly(pairwise.constraints = pairwise,
                                                     ncores = 6)
Nsamps <- 10
system.time(samples <- hitandrun(constr = all.constraints,
                                 n.samples = Nsamps))


matrices <- lapply(1:nrow(samples),
                   function(i)
                     build_matrix(vec = samples[i,],
                                  ms = mprobs))

checks <- unlist(lapply(matrices, bindata::check.commonprob))
checks