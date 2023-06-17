build_index_matrix <- function(N) {
  ind.mat <- matrix(0, nrow = N, ncol = N)
  ind.mat[lower.tri(ind.mat)] <- 1:(N * (N - 1) / 2)
  ind.mat <- ind.mat + t(ind.mat)
  diag(ind.mat) <- -1
  return(ind.mat)
}


get_pairwise_constraints <- function(ps, ncores = 6) {
  N <- length(ps)
  IDM <- build_index_matrix(N)
  the.grid <- subset(expand.grid(1:N, 1:N), Var1 < Var2)
  K <- nrow(the.grid)
  res <- parallel::mclapply(1:K,
                            function(k) {
                              ii <- the.grid[k, 1]
                              jj <- the.grid[k, 2]
                              index <- IDM[ii, jj]
                              l <- max(c(0, ps[ii] + ps[jj] - 1))
                              u <- min(c(ps[ii], ps[jj]))
                              data.frame(
                                i = ii,
                                j = jj,
                                lb = l,
                                ub = u,
                                par_index = index,
                                par_name = paste0("theta_", index)
                              )
                            },
                            mc.cores = ncores)
  out <- do.call(rbind, res)
  return(list(ind_mat = IDM,
              pairwise_ctr = out))
}

get_triplet_constraints <- function(ps) {
  P <- length(ps)
  ans.list <- vector(choose(P, 3), mode = "list")
  counter <- 1
  for (i in 1:(P - 2)) {
    for (j in (i + 1):(P - 1)) {
      for (k in (j +  1):P) {
        ans.list[[counter]] <- data.frame(
          i = i,
          j = j,
          k = k,
          l1 = ps[i] + ps[j] + ps[k] - 1,
          l2 = max(ps[i] + ps[j] - 1, 0) +
            max(ps[i] + ps[k] - 1, 0) +
            max(ps[j] + ps[k] - 1, 0),
          u = min(ps[i], ps[j]) +
            min(ps[i], ps[k]) +
            min(ps[j], ps[k])
        )
        counter <- counter + 1
      }
    }
  }
  ans <-
    tibble::tibble(do.call(rbind, ans.list))
  return(ans)
}


make_constraints_har <- function(pairwise.constraints,
                                 triplet.constraints,
                                 ncores = 6){
  
  ### Pairwise stuff
  
  prw.tab <- pairwise.constraints$pairwise_ctr
  ind.mat <- pairwise.constraints$ind_mat
  M <- nrow(prw.tab)
  
  l.prw <- parallel::mclapply(1:M, function(t){
    pp <- prw.tab[t, ]$par_index
    lowerBoundConstraint(n = M, i = pp, x = prw.tab[t, ]$lb)
  }, mc.cores = ncores)
  
  u.prw <- parallel::mclapply(1:M, function(t){
    pp <- prw.tab[t, ]$par_index
    upperBoundConstraint(n = M, i = pp, x = prw.tab[t, ]$ub)
  }, mc.cores = ncores)
  
  ## Triplet stuff
  
  trpl <- triplet.constraints
  
  Q <- nrow(trpl)
  
  l.trpl <- parallel::mclapply(1:Q, function(t){
    the.vec <- rep(0, M)
    pos <- c(
      ind.mat[trpl[t, ]$i, trpl[t, ]$j],
      ind.mat[trpl[t, ]$i, trpl[t, ]$k],
      ind.mat[trpl[t, ]$j, trpl[t, ]$k]
    )
    the.vec[pos] <- -1
    Lconstr <- list(constr = the.vec,
                    rhs = -1 * max(c(trpl[t, ]$l1, trpl[t, ]$l2)),
                    dir = "<=")
    return(Lconstr)
  }, mc.cores = ncores)
  
  u.trpl <- parallel::mclapply(1:Q, function(t){
    the.vec <- rep(0, M)
    pos <- c(
      ind.mat[trpl[t, ]$i, trpl[t, ]$j],
      ind.mat[trpl[t, ]$i, trpl[t, ]$k],
      ind.mat[trpl[t, ]$j, trpl[t, ]$k]
    )
    the.vec[pos] <- 1
    Uconstr <- list(constr = the.vec,
                    rhs = trpl[t, ]$u,
                    dir = "<=")
    return(Uconstr)
  }, mc.cores = ncores)
  
  the.ctrs <- Reduce(function(x, y) mergeConstraints(x, y),
                     c(l.prw, u.prw, l.trpl, u.trpl))
  return(the.ctrs)
}

vec2symMat <- function(x, diag = TRUE, byrow = FALSE) {
  # stolen from https://github.com/mikewlcheung/metasem/blob/17542243fc3f69c1b6d3c53ec68c00f4f8bbb81f/R/vec2symMat.R
  m <- length(x)
  d <- if (diag)
    1
  else-1
  n <- floor((sqrt(1 + 8 * m) - d) / 2)
  if (m != n * (n + d) / 2)
    stop("Cannot make a square matrix as the length of \"x\" is incorrect.")
  mat <- diag(n)
  
  ## Row major
  if (byrow) {
    mat[upper.tri(mat, diag = diag)] <- x
    index <- lower.tri(mat)
    mat[index] <- t(mat)[index]
  } else {
    ## Column major: default behavior
    mat[lower.tri(mat, diag = diag)] <- x
    # Just mirroring the matrix, exclude the diagonals
    ## mat[upper.tri(mat, diag=FALSE)] <- mat[lower.tri(mat, diag=FALSE)]
    ## Corrected a bug
    index <- upper.tri(mat)
    mat[index] <- t(mat)[index]
  }
  mat
}

build_matrix <- function(vec, ms){
  mat <- vec2symMat(x = vec, diag = FALSE, byrow = FALSE)
  diag(mat) <- ms
  return(mat)
}