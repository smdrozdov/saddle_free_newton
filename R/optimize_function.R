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

VectorIsZero <- function(v){
  if (sum(v ^ 2) < 0.00001) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

TestAll <- function(){
  point.container <- pointContainer(p = c(0, 0),
                                    L = function(v) {v[1] + 2 * v[2]},
                                    input.dimension = 2,
                                    epsilon = 0.00001,
                                    EdgeDistance = function(v) {1})
  assert_that(VectorIsZero(NumericGradient(point.container) - c(1, 2)))

  point.container <- pointContainer(p = c(0, 0),
                                    L = function(v) {v[1] ^ 2 - 2 * v[2] ^ 2},
                                    input.dimension = 2,
                                    epsilon = 0.00001,
                                    EdgeDistance = function(v) {1})
  assert_that(VectorIsZero(NumericHessian(point.container) - diag(c(2, -4))))

  point.container <- pointContainer(p = c(pi / 3, pi / 6),
                                    L = function(v) {sin(v[1]) * cos(v[2])},
                                    input.dimension = 2,
                                    epsilon = 0.00001,
                                    EdgeDistance = function(v) {1})
  assert_that(VectorIsZero(NumericHessian(point.container) - c(c(-0.75, -0.25), c(-0.25, -0.75))))

  point.container <- pointContainer(p = c(1, 0),
                                    L = function(v) {v[1] ^ 2 - v[1] * v[2] + 2 * v[2] ^ 2},
                                    input.dimension = 2,
                                    epsilon = 0.00001,
                                    EdgeDistance = function(v) {1})
  assert_that(VectorIsZero(MultiplyHessianVector(point.container, c(1, -1)) - c(3, -5)))

  result.gradient.descent <- OptimizeFunction(function(v){(v[1] + 0.1) ^ 2 + v[2] ^ 2 + (v[3] - 0.1) ^ 2},
                            500,
                            3,
                            0.00001,
                            0.000001,
                            0.1,
                            function(v) {1},
                            "GD")
  assert_that(VectorIsZero(result.gradient.descent - c(-0.1, 0, 0.1)))

  result.saddle.free.newton <- OptimizeFunction(function(v){(v[1] + 0.1) ^ 2 + v[2] ^ 2 + (v[3] - 0.1) ^ 2},
                                              500,
                                              3,
                                              0.00001,
                                              0.000001,
                                              0.1,
                                              function(v) {1},
                                              "SFN")
  assert_that(VectorIsZero(result.saddle.free.newton - c(-0.1, 0, 0.1)))

  x <- c(0.1, 0.2, 0.4, 0.5)
  y <- c(0.3, 0.4, 0.5, 0.6)
  NonConvexFunction <- function(v){
    a <- v[1]
    b <- v[2]
    c <- v[3]
    d <- v[4]
    return(sum((exp(- a * x) - exp(- b * x) - y) ^ 2) +
           sum((exp(- b * x) - exp(- c * x) - y) ^ 2) +
           sum((exp(- c * x) - exp(- d * x) - y) ^ 2) + a)
  }
  result.saddle.free.newton <- OptimizeFunction(NonConvexFunction,
                                              500,
                                              4,
                                              0.00001,
                                              0.000001,
                                              0.1,
                                              function(v) {1},
                                              "SFN")
  result.gradient.descent <- OptimizeFunction(NonConvexFunction,
                                             5000,
                                             4,
                                             0.00001,
                                             0.000001,
                                             0.1,
                                             function(v) {1},
                                             "GD")
  assert_that(NonConvexFunction(result.saddle.free.newton) < NonConvexFunction(result.gradient.descent))

  point.container <- pointContainer(p = c(1,1,1),
                                    L = function(v){return ((v[1] + 0.1) ^ 2 + v[2] ^ 2 - (v[3] - 0.1) ^ 2)},
                                    input.dimension = 3,
                                    epsilon = 0.0001,
                                    EdgeDistance = function(v) {1})
  assert_that(abs(DirectionCurvature(point.container, c(1.0, 0.0, 0.0)) - 2.0) < 0.000001)
  assert_that(abs(DirectionCurvature(point.container, c(1.0, 1.0, 0.0)) - 2.0) < 0.000001)
  assert_that(abs(DirectionCurvature(point.container, c(0.0, 0.0, 1.0)) + 2.0) < 0.000001)
  assert_that(abs(DirectionCurvature(point.container, c(1.0, 1.0, 1.0)) - 2.0 /3.0) < 0.000001)


  resASFN <- OptimizeFunction(function(v) { ((v[1] - 0.1) ^ 2 + 2 * (v[10] + 1.0) ^ 2 + 3 * (v[15] - 1.0) ^ 2) - sin(v[1] * v[15]) ^ 2},
                              50,
                              15,
                              0.00001,
                              0.000001,
                              0.3,
                              function(v) {1},
                              "ASFN",
                              2)
  print(resASFN)
}


TestKrylov <- function(){
  # TODO(smdrozdov): Finalize Krylov test.
  point.container <- pointContainer(p = c(1,1,1,1,1,1,1,1,1,1),
                                    L = function(v){return (v[1] ^ 2 - 2 * v[6] ^ 2 + v[10] ^ 2)},
                                    input.dimension = 10,
                                    epsilon = 0.0001,
                                    EdgeDistance = function(v) {1})

  KS <- KrylovSubspace(point.container, 4)$subspace
  #print(KS)
  for (i in 1:4){
    #print(DirectionCurvature(point.container, KS[i,]))
  }
  #print(DirectionCurvature(point.container, c(1,0,0,0,0,0,0,0,0,0)))
}

