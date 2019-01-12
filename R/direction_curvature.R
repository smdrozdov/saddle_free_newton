DirectionCurvature <- function(point.container,
                               direction){
  # Computes curvature in a given direction.
  # Args:
  #   point.container: point, function etc.
  #   direction: curvature direction.

  epsilon.shift <- min(point.container$epsilon.shift.input, point.container$EdgeDistance(point.container$p) / exp(1))
  direction.normalised <- direction / sqrt(sum(direction * direction))

  p.shift.positive <- point.container$p + epsilon.shift * direction.normalised
  p.shift.negative <- point.container$p - epsilon.shift * direction.normalised

  curvature <- (point.container$L(p.shift.positive)
                  - 2 * point.container$L(point.container$p)
                  + point.container$L(p.shift.negative)) / (4 * epsilon.shift ^ 2)
  return(curvature)
}
