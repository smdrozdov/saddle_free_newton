NumericHessian <- function(p,
                           L,
                           input.dimension,
                           epsilon.shift.input,
                           EdgeDistance){
  # Computes numeric hessian.
  # Args:
  #   p: initial point.
  #   L: function to be optimized, takes vector as input.
  #   input.dimension: dimension of L input.
  #   epsilon.shift.input: size of step.
  #   EdgeDistance: distance from given point to the domain boundary.

  epsilon.shift <- min(epsilon.shift.input, EdgeDistance(p) / exp(1))
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
  return(hessian)
}