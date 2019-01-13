NumericGradient <- function(point.container){
  # Computes numeric gradient.
  # Args:
  #   point.container: point and function.

  epsilon.shift <- min(point.container$epsilon, point.container$EdgeDistance(point.container$p) / exp(1))
  gradient <- vector(length = point.container$input.dimension)
  for (i in 1:point.container$input.dimension){
    p.increase.i <- point.container$p
    p.increase.i[i] <- point.container$p[i] + epsilon.shift

    p.decrease.i <- point.container$p
    p.decrease.i[i] <- point.container$p[i] - epsilon.shift
    gradient[i] <- (point.container$L(p.increase.i) - point.container$L(p.decrease.i)) / (2 * epsilon.shift)
  }
  return(gradient)
}
