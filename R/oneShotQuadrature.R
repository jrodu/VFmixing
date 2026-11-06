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

  # Weighted path-velocity design: rows (node, dim), cols (path region, dim, cubic coefficient)
  Psi = matrix(0, n_t*n_dim, 4*n_dim*Rp)
  dB = cbind(3*path$s^2, 2*path$s, 1, 0)
  for(d in seq_len(n_dim)){
    rows = (seq_len(n_t) - 1)*n_dim + d
    for(k in 1:4){
      cols = (path$idx - 1)*4*n_dim + (d - 1)*4 + k
      Psi[cbind(rows, cols)] = sw * dB[, k]
    }
  }

  list(lower = lower, upper = upper, t = t, w = w, vel = path$vel, F = Fa,
       G = crossprod(A),          # this region's block of Sigma_Traj
       S = crossprod(A, Psi))     # this region's rows of Sigma_PathTraj
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

# Transition matrix for the trajectory free parameters: the same as
# constructUnboundedCubicSplineTransitionMatrix_Traj(), built with kronecker().
trajTransitionMatrix = function(time_splits, n_models, V_Smooth){
  Tp = constructUnboundedCubicSplineTransitionMatrix_Path(time_splits, V_Smooth)
  I = diag(n_models)
  do.call(rbind, lapply(seq_len(length(time_splits) - 1), function(i){
    kronecker(I, Tp[1:4 + (i - 1)*4, , drop = FALSE])
  }))
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
