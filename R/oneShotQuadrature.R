# Internal numerical helpers for the OneShot trajectory method.
#
# OneShot needs, for every trajectory region, the integrals
#   Sigma_Traj:     int (F w-basis)' (F w-basis) dt
#   Sigma_PathTraj: int (F w-basis)' (path velocity basis) dt
# and the integrated squared velocity error. All of these share the same
# integrand ingredients: the path position and velocity, and the base vector
# fields evaluated along the path. Instead of integrating every matrix entry
# separately, these helpers evaluate the ingredients once at Gauss-Legendre
# nodes on each piece of the region (split at the path nodes, where the path
# changes cubic) and form every integral with matrix products.

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
gaussLegendre = function(k){
  if(k == 1){
    return(list(x = 0, w = 2))
  }
  i = seq_len(k - 1)
  b = i / sqrt(4*i^2 - 1)
  J = matrix(0, k, k)
  J[cbind(i, i + 1)] = b
  J[cbind(i + 1, i)] = b
  e = eigen(J, symmetric = TRUE)
  ord = order(e$values)
  list(x = e$values[ord], w = 2*e$vectors[1, ord]^2)
}

# Path region, local time, position and velocity at many times at once.
# The region rule matches getEstPathPosition(): a time on a node belongs to
# the region that ends there.
evaluatePathAt = function(t, EstPath){
  P = EstPath$TimeSplits
  R = length(P) - 1
  idx = pmin(pmax(findInterval(t, P, left.open = TRUE), 1), R)
  s = t - P[idx]
  n_dim = nrow(EstPath$Path[[1]])
  pos = vel = matrix(0, length(t), n_dim)
  for(r in unique(idx)){
    sel = idx == r
    ss = s[sel]
    C = EstPath$Path[[r]]
    pos[sel, ] = cbind(ss^3, ss^2, ss, 1) %*% t(C)
    vel[sel, ] = cbind(3*ss^2, 2*ss, 1, 0) %*% t(C)
  }
  list(idx = idx, s = s, pos = pos, vel = vel)
}

# Everything OneShot needs about one trajectory region [lower, upper] that does
# not depend on the trajectory coefficients. Computed once when the region is
# created and reused for every later candidate.
oneShotRegionCache = function(lower, upper, EstPath, baseVectorFields, n_dim, n_models, gl){

  P = EstPath$TimeSplits
  Rp = length(P) - 1

  # Quadrature nodes on each piece between path nodes
  brks = sort(unique(c(lower, upper, P[P > lower & P < upper])))
  n_pieces = length(brks) - 1
  half = diff(brks)/2
  mid = (brks[-1] + brks[-length(brks)])/2
  t = c(outer(gl$x, half) + rep(mid, each = length(gl$x)))
  w = c(outer(gl$w, half))
  n_t = length(t)
  sw = sqrt(w)

  path = evaluatePathAt(t, EstPath)

  Fa = array(0, c(n_dim, n_models, n_t))
  for(i in seq_len(n_t)){
    Fa[, , i] = baseVectorFields(t[i], path$pos[i, ])
  }

  # Weighted trajectory design: rows (node, dim), cols (model, cubic coefficient)
  u = t - lower
  B = cbind(u^3, u^2, u, 1)
  A = matrix(0, n_t*n_dim, 4*n_models)
  for(d in seq_len(n_dim)){
    rows = (seq_len(n_t) - 1)*n_dim + d
    for(m in seq_len(n_models)){
      A[rows, (m - 1)*4 + 1:4] = sw * Fa[d, m, ] * B
    }
  }

  # This region's rows of Sigma_PathTraj: S = A' Psi, with Psi the weighted
  # path-velocity design (rows (node, dim), cols (path region, dim, cubic
  # coefficient)). A node lies in one path region, so each row of Psi has only
  # 4 non-zeros; rather than build Psi (dense, its size grows with the square of
  # the number of path regions), each path region's columns of S are summed
  # directly from the nodes inside it. The same sums, without the zeros.
  #
  # Only the path regions this region overlaps (a contiguous range) have
  # non-zero columns, so only those are stored: S holds them, and S_cols gives
  # their columns among all n_path_cols path coefficients (see pathTrajFull()).
  touched = range(path$idx)
  offset = (touched[1] - 1)*4*n_dim
  S_cols = offset + seq_len((touched[2] - touched[1] + 1)*4*n_dim)
  S = matrix(0, 4*n_models, length(S_cols))
  dB = cbind(3*path$s^2, 2*path$s, 1)   # the fourth basis function's derivative is 0
  for(d in seq_len(n_dim)){
    rows = (seq_len(n_t) - 1)*n_dim + d
    A_d = A[rows, , drop = FALSE]
    for(k in 1:3){
      sums = rowsum(A_d * (sw * dB[, k]), path$idx)   # one row per path region touched
      cols = (as.integer(rownames(sums)) - 1)*4*n_dim + (d - 1)*4 + k - offset
      S[, cols] = t(sums)
    }
  }

  list(lower = lower, upper = upper, t = t, w = w, vel = path$vel, F = Fa,
       G = crossprod(A),          # this region's block of Sigma_Traj
       S = S, S_cols = S_cols, n_path_cols = 4*n_dim*Rp)
}

# A region's rows of Sigma_PathTraj across all the path coefficients (zeros
# outside the path regions it overlaps).
pathTrajFull = function(cache){
  out = matrix(0, nrow(cache$S), cache$n_path_cols)
  out[, cache$S_cols] = cache$S
  out
}

# Sigma_PathTraj as a sparse matrix: every region's rows, stacked.
pathTrajSparse = function(RegionCaches){
  k = nrow(RegionCaches[[1]]$S)
  i = unlist(lapply(seq_along(RegionCaches), function(r){rep((r - 1)*k + seq_len(k), length(RegionCaches[[r]]$S_cols))}))
  j = unlist(lapply(RegionCaches, function(cache){rep(cache$S_cols, each = k)}))
  x = unlist(lapply(RegionCaches, function(cache){c(cache$S)}))
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(k*length(RegionCaches), RegionCaches[[1]]$n_path_cols))
}

# Integrated squared velocity error over one region for given trajectory coefficients.
oneShotRegionError = function(cache, TrajCoef, n_dim, n_models){
  u = cache$t - cache$lower
  wv = TrajCoef %*% t(cbind(u^3, u^2, u, 1))          # n_models x n_t
  est = matrix(0, n_dim, length(u))
  for(m in seq_len(n_models)){
    est = est + matrix(cache$F[, m, ], nrow = n_dim) * rep(wv[m, ], each = n_dim)
  }
  sum(cache$w * colSums((t(cache$vel) - est)^2))
}

# Transition matrix for the trajectory free parameters.
trajTransitionMatrix = function(time_splits, n_models, V_Smooth){
  constructUnboundedCubicSplineTransitionMatrix_Traj(time_splits, n_models, V_Smooth)
}

# Dense block-diagonal matrix from a list of square blocks.
blockDiagDense = function(blocks){
  n = vapply(blocks, nrow, integer(1))
  off = c(0, cumsum(n))
  out = matrix(0, off[length(off)], off[length(off)])
  for(i in seq_along(blocks)){
    out[off[i] + seq_len(n[i]), off[i] + seq_len(n[i])] = blocks[[i]]
  }
  out
}

# pinv(A) %*% B for a symmetric positive semi-definite A. Uses a Cholesky solve
# when A is comfortably well conditioned, which gives the same answer as the
# pseudo-inverse; otherwise falls back to pracma::pinv() exactly as before.
psdSolve = function(A, B){
  ch = tryCatch(chol(A), error = function(e){NULL})
  if(!is.null(ch)){
    dg = diag(ch)
    # diag(chol)^2 tracks the eigenvalue spread; require condition number below ~1e12
    if(min(dg) > 1e-6*max(dg)){
      return(backsolve(ch, forwardsolve(t(ch), B)))
    }
  }
  pracma::pinv(A) %*% B
}

# blockdiag(blocks) %*% M without forming the block-diagonal matrix.
blockDiagTimes = function(blocks, M){
  n = vapply(blocks, nrow, integer(1))
  off = c(0, cumsum(n))
  do.call(rbind, lapply(seq_along(blocks), function(i){blocks[[i]] %*% M[off[i] + seq_len(n[i]), , drop = FALSE]}))
}

# M %*% blockdiag(blocks) without forming the block-diagonal matrix.
timesBlockDiag = function(M, blocks){
  n = vapply(blocks, nrow, integer(1))
  off = c(0, cumsum(n))
  do.call(cbind, lapply(seq_along(blocks), function(i){M[, off[i] + seq_len(n[i]), drop = FALSE] %*% blocks[[i]]}))
}

# The OneShot trajectory coefficients for C1 trajectories, without the
# covariance, from a banded solve. The least-squares fit depends only on the
# spline space, not on its basis, so this uses the local cubic Hermite basis:
# each trajectory's value and slope at every node. Each region then touches
# only its two end nodes, the normal matrix is banded, and a sparse Cholesky
# solve costs about linear time in the number of regions (the cumulative basis
# of trajTransitionMatrix() makes it dense). The coefficients are returned in
# the usual per-region (dt^3, dt^2, dt, 1) form, the same layout as
# runEstTrajStep_OneShot() uses. NULL if the system is singular or badly
# conditioned, so the caller can fall back to the dense solve.
oneShotCoefficientsC1 = function(RegionCaches, TimeSplits, path_coef, n_models){

  n_regions = length(TimeSplits) - 1
  m2 = 2*n_models          # unknowns per node: a value and a slope per trajectory
  k = 4*n_models           # unknowns per region (its two end nodes), and its coefficients

  # Hermite map from (y0, s0, y1, s1) to the cubic's (dt^3, dt^2, dt, 1) coefficients
  hermite = function(h){
    rbind(c(2/h^3, 1/h^2, -2/h^3, 1/h^2),
          c(-3/h^2, -2/h, 3/h^2, -1/h),
          c(0, 1, 0, 0),
          c(1, 0, 0, 0))
  }

  # Unknowns are ordered by node, then trajectory, then (value, slope)
  local = function(m){c(2*(m - 1) + 1:2, m2 + 2*(m - 1) + 1:2)}

  blocks = lapply(seq_len(n_regions), function(r){
    H = hermite(TimeSplits[r + 1] - TimeSplits[r])
    L = matrix(0, k, k)
    for(m in seq_len(n_models)){
      L[(m - 1)*4 + 1:4, local(m)] = H
    }
    list(L = L,
         K = crossprod(L, RegionCaches[[r]]$G %*% L),
         b = crossprod(L, RegionCaches[[r]]$S %*% path_coef[RegionCaches[[r]]$S_cols]))
  })

  n_unknowns = m2*(n_regions + 1)
  offsets = (seq_len(n_regions) - 1)*m2
  idx = expand.grid(i = seq_len(k), j = seq_len(k))

  normal = Matrix::sparseMatrix(i = unlist(lapply(offsets, function(o){o + idx$i})),
                                j = unlist(lapply(offsets, function(o){o + idx$j})),
                                x = unlist(lapply(blocks, function(bl){c(bl$K)})),
                                dims = c(n_unknowns, n_unknowns))
  normal = Matrix::forceSymmetric(normal)

  rhs = numeric(n_unknowns)
  for(r in seq_len(n_regions)){
    rows = offsets[r] + seq_len(k)
    rhs[rows] = rhs[rows] + blocks[[r]]$b
  }

  # Same conditioning rule as psdSolve(): otherwise let the dense solve handle it
  ch = tryCatch(Matrix::chol(normal), error = function(e){NULL}, warning = function(w){NULL})
  if(is.null(ch)){
    return(NULL)
  }
  dg = Matrix::diag(ch)
  if(!(min(dg) > 1e-6*max(dg))){
    return(NULL)
  }

  theta = as.numeric(Matrix::solve(ch, Matrix::solve(Matrix::t(ch), rhs)))

  unlist(lapply(seq_len(n_regions), function(r){
    coef = blocks[[r]]$L %*% theta[offsets[r] + seq_len(k)]
    c(coef)
  }))

}

# The diagonal block of the trajectory coefficient covariance for each region,
# M Cov(path) M', without forming the whole matrix. The path covariance is
# block-diagonal, one block per dimension, so each dimension's columns of M are
# used with its own block, and only the regions' diagonal blocks are formed.
trajCovBlocks = function(M, path_cov, n_regions, k){

  dims_ok = is.list(path_cov) && sum(vapply(path_cov, nrow, integer(1))) == ncol(M)

  if(!dims_ok){
    full = M %*% tcrossprod(as.matrix(Matrix::bdiag(path_cov)), M)
    return(lapply(seq_len(n_regions), function(x){full[(x - 1)*k + seq_len(k), (x - 1)*k + seq_len(k)]}))
  }

  blocks = lapply(seq_len(n_regions), function(x){matrix(0, k, k)})
  start = 0
  for(C in path_cov){
    cols = start + seq_len(nrow(C))
    M_d = M[, cols, drop = FALSE]
    W_d = M_d %*% C
    for(x in seq_len(n_regions)){
      rows = (x - 1)*k + seq_len(k)
      blocks[[x]] = blocks[[x]] + tcrossprod(W_d[rows, , drop = FALSE], M_d[rows, , drop = FALSE])
    }
    start = start + nrow(C)
  }
  blocks

}
