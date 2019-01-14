# Unittests.
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

  f <- function(v) { sum((v + 0.1) ^ 4) + v[3] ^ 6 + v[4] ^ 6 + v[5] ^ 6 - sin(v[7] * v[8]) - cos(v[1] + v[3] + v[5] + v[12] + pi / 3)}
  resGD <- OptimizeFunction(f,
                            3000,
                            15,
                            0.00001,
                            0.0000000001,
                            0.3,
                            function(v) {1},
                            "GD",
                            2)
  resASFN <- OptimizeFunction(f,
                              300,
                              15,
                              0.00001,
                              0.0000000001,
                              0.7,
                              function(v) {1},
                              "ASFN",
                              2)
  print(resASFN)
  print(resGD)
  print(f(resASFN))
  print(f(resGD))
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
