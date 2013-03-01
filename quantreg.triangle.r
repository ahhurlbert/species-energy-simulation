# Function for plotting quantile regression triangle between two quantile regression lines,
# with the assumption that the variation in the y variable is narrow at low values of x
# and higher for high values of x. Shades the triangle.

quantreg.triangle = function(x, y, upper.quant, lower.quant) {
  require(quantreg)
  qr.U = rq(y~x, upper.quant)
  qr.L = rq(y~x, lower.quant)
  slope.U = qr.U$coefficients[2]
  slope.L = qr.L$coefficients[2]
  int.U = qr.U$coefficients[1]
  int.L = qr.L$coefficients[1]
  #Find the intersection point of the two lines
  x.ints = (int.L - int.U)/(slope.U - slope.L)
  y.ints = (slope.U*int.L - slope.L*int.U)/(slope.U - slope.L)
  y.U = slope.U*max(x) + int.U
  y.L = slope.L*max(x) + int.L
  polygon(c(x.ints, max(x), max(x)), c(y.ints, y.U, y.L), col = rgb(.1, .1, .1, .1))
  
}

