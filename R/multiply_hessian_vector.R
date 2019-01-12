MultiplyHessianVector <- function(point.container,
                                  v){
  # Computes product of hessian and arbitrary vector, Hv, in O(d)-time, where d is domain dimension.
  # See http://www.bcl.hamilton.ie/~barak/papers/nc-hessian.pdf.
  # Args:
  #   point.container: point and function.
  #   v: vector to multiply.
  epsilon.shift <- min(point.container$epsilon.shift.input, point.container$EdgeDistance(point.container$p) / exp(1))
  gradient <- NumericGradient(point.container)
  delta <- epsilon.shift * v
  point.container$p <- as.vector(point.container$p + delta)
  gradient.shift <- NumericGradient(point.container)
  point.container$p <- as.vector(point.container$p - delta)

  return ((gradient.shift - gradient) / epsilon.shift)
}
