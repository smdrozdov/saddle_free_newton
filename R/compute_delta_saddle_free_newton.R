ComputeDeltaSaddleFreeNewton <- function(p,
                                         L,
                                         input.dimension,
                                         epsilon.shift.input,
                                         EdgeDistance){
  # Computes direction of maximal descent with saddle-free Newton method.
  # Args:
  #   p: initial point.
  #   L: function to be optimized, takes vector as input.
  #   input.dimension: dimension of L input.
  #   epsilon.shift.input: size of step.
  #   EdgeDistance: distance from given point to the domain boundary.

  epsilon.shift <- min(epsilon.shift.input, EdgeDistance(p) / exp(1))
  gradient <- NumericGradient(p, L, input.dimension, epsilon.shift.input, EdgeDistance)
  hessian <- NumericHessian(p, L, input.dimension, epsilon.shift.input, EdgeDistance)

  # Compute hessian absolute value, |H| in Ganguli's notation.
  ev <- eigen(hessian)
  vectors <- ev$vectors
  values <- ev$values
  hessian.pos <- vectors %*% diag(abs(values)) %*% t(vectors)

  delta <- - gradient %*% solve(hessian.pos)
  return(delta)
}
