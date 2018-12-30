NumericGradient <- function(p,
                            L,
                            input.dimension,
                            epsilon.shift.input,
                            EdgeDistance){
  # Computes numeric gradient.
  # Args:
  #   p: initial point.
  #   L: function to be optimized, takes vector as input.
  #   input.dimension: dimension of L input.
  #   epsilon.shift.input: size of step.
  #   EdgeDistance: distance from given point to the domain boundary.

  epsilon.shift <- min(epsilon.shift.input, EdgeDistance(p) / exp(1))
  gradient <- vector(length = input.dimension)
  for (i in 1:input.dimension){
    p.increase.i <- p
    p.increase.i[i] <- p[i] + epsilon.shift

    p.decrease.i <- p
    p.decrease.i[i] <- p[i] - epsilon.shift
    gradient[i] <- (L(p.increase.i) - L(p.decrease.i)) / (2 * epsilon.shift)
  }
  return(gradient)
}