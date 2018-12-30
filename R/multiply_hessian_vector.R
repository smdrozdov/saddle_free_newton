MultiplyHessianVector <- function(p,
                                  L,
                                  input.dimension,
                                  epsilon.shift.input,
                                  EdgeDistance,
                                  v){
  # Computes product of hessian and arbitrary vector, Hv, in O(d)-time, where d is domain dimension.
  # See http://www.bcl.hamilton.ie/~barak/papers/nc-hessian.pdf.
  # Args:
  #   p: initial point.
  #   L: function to be optimized, takes vector as input.
  #   input.dimension: dimension of L input.
  #   epsilon.shift.input: size of step.
  #   EdgeDistance: distance from given point to the domain boundary.
  #   v: vector to multiply.
  epsilon.shift <- min(epsilon.shift.input, EdgeDistance(p) / exp(1))

  gradient <- NumericGradient(p, L, input.dimension, epsilon.shift.input, EdgeDistance)
  gradient.shift <- NumericGradient(p + epsilon.shift * v, L, input.dimension, epsilon.shift.input, EdgeDistance)

  return ((gradient.shift - gradient) / epsilon.shift)
}
