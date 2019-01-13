KrylovSubspace <- function(point.container,
                           k){
  # Computes Krylov subspace of dimension k.
  # Args:
  #   point.container: point and function.
  #   k: Krylov subspace dimension.

  subspace <- vector(length = point.container$input.dimension)
  subspace.multiplied.by.hessian <- vector(length = point.container$input.dimension)

  g <- - NumericGradient(point.container)
  v <- g / sqrt(sum(g * g))

  subspace <- rbind(subspace, v)
  subspace.multiplied.by.hessian <- rbind(subspace.multiplied.by.hessian, MultiplyHessianVector(point.container, v))
  beta <- 0
  for (i in 2:k){
    w <- MultiplyHessianVector(point.container, v)
    alpha <- sum(w * v)
    w <- w - alpha * subspace[i,] - beta * subspace[i - 1,]
    beta <- sqrt(sum(w * w))
    v <- w / sqrt(sum(w * w))
    subspace <- rbind(subspace, v)
    subspace.multiplied.by.hessian <- rbind(subspace.multiplied.by.hessian, MultiplyHessianVector(point.container, v))
  }
  res <- KrylovContainer(subspace = subspace[-1,],
                         subspace.multiplied.by.hessian = subspace.multiplied.by.hessian[-1,])
  return(res)
}
