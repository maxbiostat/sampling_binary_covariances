# Sampling random covariance matrices for multivariate binary distributions

See [notes](https://github.com/maxbiostat/sampling_binary_covariances/blob/main/notes/sampling_binary_covariance.pdf) for the problem description. 

Take **X** a binary vector of dimension P, with some distribution `F(x ; theta)`. We would like to sample from the space of cross-moment (joint probability) matrices uniformly given only a restriction on the expectation of **X**
In short, given a P-dimensional vector of marginal probabilities **m**, we want to sample a matrix **C** (uniformly) from the space of covariance matrices compatible with **m**, Q(**m**). 
