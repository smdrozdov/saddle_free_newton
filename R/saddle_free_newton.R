# Implementation of Ganguli-Bengio method of non-convex function optimzation.
# Link to article: https://arxiv.org/abs/1406.2572
#   TODO(smdrozdov): Introduce default values.
#   TODO(smdrozdov): Lanczosh vector procedures.
# Args:
#   L: function to be optimized, takes vector as input.
#   max.steps: amount of steps, set to 500 by default.
#   input.dimension: dimension of L input.
#   epsilon.shift.input: size of step.
#   epsilon.stop: stop in this case.
#   learning.rate: 0.1 by default.
#   EdgeDistance: distance from given point to the domain boundary.
SaddleFreeNewton <- function(L,
                             max.steps,
                             input.dimension,
                             epsilon.shift.input,
                             epsilon.stop,
                             learning.rate,
                             EdgeDistance){
  shifts <- {}
  p <- vector(length = input.dimension)
  for (step in 1:max.steps){
    gradient <- vector(length = input.dimension)
    epsilon.shift <- min(epsilon.shift.input, EdgeDistance(p) / exp(1))
    for (i in 1:input.dimension){
      p.increase.i <- p
      p.increase.i[i] <- p[i] + epsilon.shift

      p.decrease.i <- p
      p.decrease.i[i] <- p[i] - epsilon.shift
      gradient[i] <- (L(p.increase.i) - L(p.decrease.i)) / (2 * epsilon.shift)
    }

    hessian <- matrix(ncol = input.dimension, nrow = input.dimension)
    # Compute numeric hessian.
    for (i in 1:input.dimension){
      for (j in 1:input.dimension){
        if (i == j){
          p.increase.i <- p
          p.increase.i[i] <- p[i] + epsilon.shift
          p.decrease.i <- p
          p.decrease.i[i] <- p[i] - epsilon.shift
          hessian[i, i] <- (L(p.increase.i) - 2 * L(p) + L(p.decrease.i)) / (epsilon.shift ^ 2)
        } else if (i > j){
          p.increase.i.j <- p
          p.increase.i.j[i] <- p[i] + epsilon.shift
          p.increase.i.j[j] <- p[j] + epsilon.shift

          p.increase.i.decrease.j <- p
          p.increase.i.decrease.j[i] <- p[i] + epsilon.shift
          p.increase.i.decrease.j[j] <- p[j] - epsilon.shift

          p.decrease.i.increase.j <- p
          p.decrease.i.increase.j[i] <- p[i] - epsilon.shift
          p.decrease.i.increase.j[j] <- p[j] + epsilon.shift

          p.decrease.i.j <- p
          p.decrease.i.j[i] <- p[i] - epsilon.shift
          p.decrease.i.j[j] <- p[j] - epsilon.shift

          hessian[i, j] <- (L(p.increase.i.j)
                            - L(p.increase.i.decrease.j)
                            - L(p.decrease.i.increase.j)
                            + L(p.decrease.i.j)) / (4 * epsilon.shift ^ 2)
        }
      }
    }

    # Add missing values below diagonal.
    for (i in 1:input.dimension){
      for (j in 1:input.dimension){
        if (i < j) {
          hessian[i, j] <- hessian[j, i]
        }
      }
    }

    # Compute hessian absolute value, |H| in Ganguli's notation.
    ev <- eigen(hessian)
    vectors <- ev$vectors
    values <- ev$values
    hessian.pos <- vectors %*% diag(abs(values)) %*% t(vectors)

    delta <- - gradient %*% solve(hessian.pos) * learning.rate
    p <- p + delta

    shifts <- c(sqrt(delta %*% t(delta)), shifts)
    s <- sum(shifts[1:10])
    if (!is.na(s) && s < epsilon.stop) {
      break
    }
  }
  return(p)
}


# Mock test of Sadddle Free Newton.
Main <- function(){
  # Lambda. Can be anonymous.
  DistanceToPoint <- function(v){
    return ((v[1] + 0.1) ^ 2 + v[2] ^ 2 + (v[3] - 0.1) ^ 2)
  }
  CircleDistance <- function(v){
    return (1 - sqrt(v[1] ^ 2 + v[2] ^ 2 + v[3] ^ 2))
  }
  one <- function(v){ return (1)}
  res <- SaddleFreeNewton(DistanceToPoint, 500, 3, 0.00001, 0.000001, 0.1, one)
  cat(res)
}


