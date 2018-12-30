# Implementation of Ganguli-Bengio method of non-convex function optimzation.
# Link to article: https://arxiv.org/abs/1406.2572
#   TODO(smdrozdov): Introduce default values.
#   TODO(smdrozdov): Approximate Saddle-Free Newton:
#                    1. Project into Krylov subspace in optimize.funciton.
#   TODO(smdrozdov): Move tests to new file.
#   TODO(smdrozdov): Test on non-convex function.

OptimizeFunction <- function(L,
                             max.steps,
                             input.dimension,
                             epsilon.shift.input,
                             epsilon.stop,
                             learning.rate,
                             EdgeDistance,
                             optimization.method){
  # Calculates global minimum of function defined in euclidian space.
  # Args:
  #   L: function to be optimized, takes vector as input.
  #   max.steps: amount of steps, set to 500 by default.
  #   input.dimension: dimension of L input.
  #   epsilon.shift.input: size of step.
  #   epsilon.stop: stop in this case.
  #   learning.rate: 0.1 by default.
  #   EdgeDistance: distance from given point to the domain boundary.

  p <- vector(length = input.dimension)
  if (optimization.method == "GD"){
    # Gradient descent.
    shifts <- {}
    for (step in 1:max.steps){
      delta <- -t(NumericGradient(p, L, input.dimension, epsilon.shift.input, EdgeDistance))
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
      delta <- ComputeDeltaSaddleFreeNewton(p, L, input.dimension, epsilon.shift.input, EdgeDistance)
      p <- p + delta * learning.rate
      shifts <- c(sqrt(sum(delta * delta)), shifts)
      s <- sum(shifts[1:10])
      if (!is.na(s) && s < epsilon.stop) {
        break
      }
    }
  } else if (optimization.method == "ASFN"){
    # Approximate Saddle-free Newton.
    #TODO(smdrozdov): Krylov subspace projection.
  }
  return(p)
}


# Mock test of Sadddle Free Newton.
TestSFN <- function(){
  # Lambda. Can be anonymous.
  DistanceToPoint <- function(v){
    return ((v[1] + 0.1) ^ 2 + v[2] ^ 2 + (v[3] - 0.1) ^ 2)
  }
  CircleDistance <- function(v){
    return (1 - sqrt(v[1] ^ 2 + v[2] ^ 2 + v[3] ^ 2))
  }
  one <- function(v){ return (1)}
  resGD <- OptimizeFunction(DistanceToPoint, 500, 3, 0.00001, 0.000001, 0.1, one, "GD")
  print(resGD)

  resSFN <- OptimizeFunction(DistanceToPoint, 500, 3, 0.00001, 0.000001, 0.1, one, "SFN")
  print(resSFN)

  Hv <- MultiplyHessianVector(c(1,1,1), DistanceToPoint, 3, 0.0001, one, c(1,1,-1))
  print(Hv)

  KS <- KrylovSubspace(c(1,1,1), DistanceToPoint, 3, 0.0001, one, c(1,1,1), 3)
  print(KS)
}


