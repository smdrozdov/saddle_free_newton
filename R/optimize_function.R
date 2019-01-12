# Implementation of Ganguli-Bengio method of non-convex function optimzation.
# Link to article: https://arxiv.org/abs/1406.2572
# Optimization methods available: GD - Gradient Descent. Here mostly for comparison.
#                                 SFN - Saddle-free Newton. Escapes saddle points, takes O(d^3),
#                                       where d is domain dimension.
#                                 ASFN - Approximate Saddle-free Newton. Escapes saddle points,
#                                        takes O(d), where d is domain dimension.
#   P0 TODO(smdrozdov): Port to Python.
#   P0 TODO(smdrozdov): Move tests to new file.
#   P1 TODO(smdrozdov): Test on non-convex function.
#   P2 TODO(smdrozdov): Approximate Saddle-Free Newton:
#                    1. Project into Krylov subspace in optimize.funciton.
#   TODO(smdrozdov): Krylov subspace unstable if gradient is exactly zero.
#   TODO(smdrozdov): Switch to Python/Tensorflow.
#   TODO(smdrozdov): Make alpha number, not matrix.
#   TODO(smdrozdov): Switch from point.container to function.container.
#   TODO(smdrozdov): Remove ComputeDeltaSaddleFreeNewton.

OptimizeFunction <- function(L,
                             max.steps,
                             input.dimension,
                             epsilon.shift.input,
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
  #   epsilon.shift.input: size of step.
  #   epsilon.stop: stop in this case.
  #   learning.rate: 0.1 by default.
  #   EdgeDistance: distance from given point to the domain boundary.
  #   optimization.mehtod: GD, SFN or ASFN.
  #   k: Krylov subspace dimension in ASFN.
  #   m: amount of subcycle iterations in ASFN.

  p <- numeric(input.dimension)
  if (optimization.method == "GD"){
    # Gradient descent.
    shifts <- {}
    for (step in 1:max.steps){
      p <- as.vector(p)
      point.container <- pointContainer(p = p,
                                        L = L,
                                        input.dimension = input.dimension,
                                        epsilon.shift.input = epsilon.shift.input,
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
                                        epsilon.shift.input = epsilon.shift.input,
                                        EdgeDistance = EdgeDistance)

      epsilon.shift <- min(point.container$epsilon.shift.input, point.container$EdgeDistance(point.container$p) / exp(1))
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
    one <- function(v){ return (1)}
    delta = vector(length = input.dimension)
    for (j in 1:k){
      delta[j] = 1
    }
    point.container <- pointContainer(p = p,
                                      L = L,
                                      input.dimension = input.dimension,
                                      epsilon.shift.input = epsilon.shift.input,
                                      EdgeDistance = EdgeDistance)
    for (step in 1:max.steps){
      p <- as.vector(p)
      point.container$p <- p
      krylov.subspace <- KrylovSubspace(point.container, delta, k)
      V <- krylov.subspace$subspace
      V.multiplied.by.hessian <- krylov.subspace$subspace.multiplied.by.hessian

      L_cap <- function(alpha){
        return(L(p + alpha %*% V))
      }

      hessian.subspace <- V.multiplied.by.hessian %*% t(V)

      # Compute hessian absolute value.
      ev <- eigen(hessian.subspace)
      vectors <- ev$vectors
      values <- ev$values
      hessian.subspace.pos <- vectors %*% diag(abs(values)) %*% t(vectors)

      id <- vector(length = k)
      for (j in 1:k){
        id[j] = 1
      }
      for (i in 1:m){
        point.container$p <- as.vector(p)
        g <- - V %*% NumericGradient(point.container)
        projection <- function(lambda){
          return(L_cap(t(g) %*% solve(hessian.subspace.pos + lambda[1] * diag(id))))
        }
        print(projection(c(0,0,0)))
        print(projection(c(0,1,0)))
        # Here complex number appear due to matrix inversion.
        lambda.min <- OptimizeFunction(projection, 10, 1, 0.00001, 0.000001, 0.1, one, "GD")
        delta <- t(g) %*% solve(hessian.subspace.pos + lambda.min[1] * diag(id)) %*% V
        p <- p + delta
      }
    }

  }
  return(p)
}


# Mock test of Sadddle Free Newton.
TestSFN <- function(){
  # Lambda. Can be anonymous.
  DistanceToPoint <- function(v){
    return ((v[1] + 0.1) ^ 2 + v[2] ^ 2 + (v[3] - 0.1) ^ 2)
  }
  DistanceToOrigin <- function(v){
    return (v[1] ^ 2 + v[2] ^ 2 + v[3] ^ 2 + v[4] ^ 2 + v[5] ^ 2)
  }
  CircleDistance <- function(v){
    return (1 - sqrt(v[1] ^ 2 + v[2] ^ 2 + v[3] ^ 2))
  }
  one <- function(v){ return (1)}
  resGD <- OptimizeFunction(DistanceToPoint, 500, 3, 0.00001, 0.000001, 0.1, one, "GD")
  print(resGD)

  resSFN <- OptimizeFunction(DistanceToPoint, 500, 3, 0.00001, 0.000001, 0.1, one, "SFN")
  print(resSFN)
}

TestDirectionCurvature <- function(){
  DistanceToPoint <- function(v){
    return ((v[1] + 0.1) ^ 2 + v[2] ^ 2 + (v[3] - 0.1) ^ 2)
  }
  one <- function(v){ return (1)}
  point.container <- pointContainer(p = c(1,1,1),
                                    L = DistanceToPoint,
                                    input.dimension = 3,
                                    epsilon.shift.input = 0.0001,
                                    EdgeDistance = one)

  direction.curvature <- DirectionCurvature(point.container, c(1.0, 1.0, 1.0))
  print(direction.curvature)
}

TestKrylov <- function(){
  DistanceToPoint <- function(v){
    return ((v[1] + 0.1) ^ 2 + v[2] ^ 2 + (v[3] - 0.1) ^ 2)
  }
  one <- function(v){ return (1)}
  point.container <- pointContainer(p = c(1,1,1),
                                    L = DistanceToPoint,
                                    input.dimension = 3,
                                    epsilon.shift.input = 0.0001,
                                    EdgeDistance = one)

  Hv <- MultiplyHessianVector(point.container, c(1,1,-1))
  print(Hv)

  KS <- KrylovSubspace(point.container, c(1,1,1), 3)$subspace
  print(KS)
}

# Mock test of Approximate Saddle-free Newton.
TestASFN <- function(){
  # Lambda. Can be anonymous.
  DistanceToPoint <- function(v){
    return ((v[1] - 1.0) ^ 2 + v[2] ^ 2 + (v[3] + 1.0) ^ 2)
  }

  one <- function(v){ return (1)}

  resSFN <- OptimizeFunction(DistanceToPoint, 3, 3, 0.00001, 0.000001, 0.1, one, "ASFN", 2, 2)
  print(resSFN)
}

