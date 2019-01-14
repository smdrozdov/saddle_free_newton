# Implementation of Ganguli-Bengio method of non-convex function optimzation.
# Link to article: https://arxiv.org/abs/1406.2572
# Optimization methods available: GD - Gradient Descent. Here mostly for comparison.
#                                 SFN - Saddle-free Newton. Escapes saddle points, takes O(d^3),
#                                       where d is domain dimension.
#                                 ASFN - Approximate Saddle-free Newton. Escapes saddle points,
#                                        takes O(d), where d is domain dimension.
#   P0 TODO(smdrozdov): Port to Python.
#   P1 TODO(smdrozdov): Test Krylov subspace with Direction curvature.

#   TODO(smdrozdov): Krylov subspace unstable if gradient is exactly zero.
#   TODO(smdrozdov): Switch from point.container to function.container.
#   TODO(smdrozdov): Add previous delta to Krylov.
#   TODO(smdrozdov): Optimize in the subspace.

OptimizeFunction <- function(L,
                             max.steps,
                             input.dimension,
                             epsilon,
                             epsilon.stop,
                             learning.rate,
                             EdgeDistance,
                             optimization.method,
                             k,
                             m){
  # Calculates global minimum of function defined in euclidian space.
  # Args:
  #   L: function to be optimized, takes vector as input.
  #   max.steps: amount of steps, set to 500 by default.
  #   input.dimension: dimension of L input.
  #   epsilon: size of step.
  #   epsilon.stop: stop in this case.
  #   learning.rate: 0.1 by default.
  #   EdgeDistance: distance from given point to the domain boundary.
  #   optimization.mehtod: GD, SFN or ASFN.
  #   k: Krylov subspace dimension in ASFN.

  p <- numeric(input.dimension)
  if (optimization.method == "GD"){
    # Gradient descent.
    shifts <- {}
    for (step in 1:max.steps){
      p <- as.vector(p)
      point.container <- pointContainer(p = p,
                                        L = L,
                                        input.dimension = input.dimension,
                                        epsilon = epsilon,
                                        EdgeDistance = EdgeDistance)
      delta <- -t(NumericGradient(point.container))
      p <- p + delta * learning.rate
      shifts <- c(sqrt(sum(delta * delta)), shifts)
      s <- sum(shifts[1:10])
      if (!is.na(s) && s < epsilon.stop) {
        break
      }
    }
  } else if (optimization.method == "SFN"){
    # Saddle-free Newton.
    shifts <- {}
    for (step in 1:max.steps){
      p <- as.vector(p)
      point.container <- pointContainer(p = p,
                                        L = L,
                                        input.dimension = input.dimension,
                                        epsilon = epsilon,
                                        EdgeDistance = EdgeDistance)

      epsilon.shift <- min(point.container$epsilon, point.container$EdgeDistance(point.container$p) / exp(1))
      gradient <- NumericGradient(point.container)
      hessian <- NumericHessian(point.container)

      # Compute hessian absolute value, |H| in Ganguli's notation.
      ev <- eigen(hessian)
      vectors <- ev$vectors
      values <- ev$values
      hessian.pos <- vectors %*% diag(abs(values)) %*% t(vectors)

      delta <- - gradient %*% solve(hessian.pos)

      p <- p + delta * learning.rate
      shifts <- c(sqrt(sum(delta * delta)), shifts)
      s <- sum(shifts[1:10])
      if (!is.na(s) && s < epsilon.stop) {
        break
      }
    }
  } else if (optimization.method == "ASFN"){
    # Approximate Saddle-free Newton.
    shifts <- {}
    one <- function(v){ return (1)}

    for (step in 1:max.steps){
      p <- as.vector(p)
      point.container <- pointContainer(p = p,
                                        L = L,
                                        input.dimension = input.dimension,
                                        epsilon = epsilon,
                                        EdgeDistance = EdgeDistance)
      gradient <- NumericGradient(point.container)

      krylov.subspace <- KrylovSubspace(point.container, k)
      V <- krylov.subspace$subspace
      V.multiplied.by.hessian <- krylov.subspace$subspace.multiplied.by.hessian
      hessian.subspace <- V.multiplied.by.hessian %*% t(V)

      # Compute hessian absolute value.
      ev <- eigen(hessian.subspace)
      vectors <- ev$vectors
      values <- ev$values
      hessian.subspace.pos <- vectors %*% diag(abs(values)) %*% t(vectors)

      delta <- - t(V %*% gradient) %*% solve(hessian.subspace.pos) %*% V

      p <- p + delta * learning.rate
      shifts <- c(sqrt(sum(delta * delta)), shifts)
      s <- sum(shifts[1:10])
      if (!is.na(s) && s < epsilon.stop) {
        break
      }
    }

  }
  return(p)
}
