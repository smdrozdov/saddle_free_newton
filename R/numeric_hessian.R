NumericHessian <- function(point.container){
  # Computes numeric hessian.
  # Args:
  #   point.container: point and function.

  epsilon.shift <- min(point.container$epsilon.shift.input, point.container$EdgeDistance(point.container$p) / exp(1))
  hessian <- matrix(ncol = point.container$input.dimension, nrow = point.container$input.dimension)
  # Compute numeric hessian.
  for (i in 1:point.container$input.dimension){
    for (j in 1:point.container$input.dimension){
      if (i == j){
        p.increase.i <- point.container$p
        p.increase.i[i] <- point.container$p[i] + epsilon.shift
        p.decrease.i <- point.container$p
        p.decrease.i[i] <- point.container$p[i] - epsilon.shift
        hessian[i, i] <- (point.container$L(p.increase.i)
                          - 2 * point.container$L(point.container$p)
                          + point.container$L(p.decrease.i)) / (epsilon.shift ^ 2)
      } else if (i > j){
        p.increase.i.j <- point.container$p
        p.increase.i.j[i] <- point.container$p[i] + epsilon.shift
        p.increase.i.j[j] <- point.container$p[j] + epsilon.shift

        p.increase.i.decrease.j <- point.container$p
        p.increase.i.decrease.j[i] <- point.container$p[i] + epsilon.shift
        p.increase.i.decrease.j[j] <- point.container$p[j] - epsilon.shift

        p.decrease.i.increase.j <- point.container$p
        p.decrease.i.increase.j[i] <- point.container$p[i] - epsilon.shift
        p.decrease.i.increase.j[j] <- point.container$p[j] + epsilon.shift

        p.decrease.i.j <- point.container$p
        p.decrease.i.j[i] <- point.container$p[i] - epsilon.shift
        p.decrease.i.j[j] <- point.container$p[j] - epsilon.shift

        hessian[i, j] <- (point.container$L(p.increase.i.j)
                          - point.container$L(p.increase.i.decrease.j)
                          - point.container$L(p.decrease.i.increase.j)
                          + point.container$L(p.decrease.i.j)) / (4 * epsilon.shift ^ 2)
      }
    }
  }

  # Add missing values below diagonal.
  for (i in 1:point.container$input.dimension){
    for (j in 1:point.container$input.dimension){
      if (i < j) {
        hessian[i, j] <- hessian[j, i]
      }
    }
  }
  return(hessian)
}
