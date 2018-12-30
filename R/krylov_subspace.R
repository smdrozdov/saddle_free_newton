KrylovSubspace <- function(p,
                           L,
                           input.dimension,
                           epsilon.shift.input,
                           EdgeDistance,
                           deltaT,
                           k){
  # Computes Krylov subspace of dimension k.
  # Args:
  #   p: initial point.
  #   L: function to be optimized, takes vector as input.
  #   input.dimension: dimension of L input.
  #   epsilon.shift.input: size of step.
  #   EdgeDistance: distance from given point to the domain boundary.
  #   k: Krylov subspace dimension.

  res <- vector(length = input.dimension)

  g <- NumericGradient(p, L, input.dimension, epsilon.shift.input, EdgeDistance)
  v <- vector(length = input.dimension)
  v <- g / sqrt(sum(g * g))
  res <- rbind(res, v)
  beta <- 0
  for (i in 2:k){
    w <- MultiplyHessianVector(p, L, input.dimension, epsilon.shift.input, EdgeDistance, v)
    if (i == k){
      w <- deltaT
    }
    alpha <- sum(w * v)
    w <- w - alpha * res[i,] - beta * res[i - 1,]
    beta <- sqrt(sum(w * w))
    v <- w / sqrt(sum(w * w))
    res <- rbind(res, v)
  }
  return (res[-1,])
}
