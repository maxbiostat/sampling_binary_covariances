source("aux.r")
library(hitandrun)
library(bindata)
library(cladeCorrelation)
###########################

ntaxa <- 4
allCs <- make_all_clades(n = ntaxa)
P <- length(allCs)
mprobs <- pclade(k = get_clade_size(allCs), n = ntaxa)

pairwise <- get_pairwise_constraints(ps = mprobs)
triplet <- get_triplet_constraints(ps = mprobs)

all.constraints <- make_constraints_har(pairwise.constraints = pairwise,
                                        triplet.constraints = triplet,
                                        ncores = 6)
Nsamps <- 10
samples <- hitandrun(constr = all.constraints,
                     n.samples = Nsamps)

matrices <- lapply(1:nrow(samples),
                   function(i) build_matrix(vec = samples[i, ],
                                            ms = mprobs))

checks <- unlist(
  lapply(matrices, bindata::check.commonprob)
)
checks