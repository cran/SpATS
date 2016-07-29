bbase <-
function(X., XL., XR., NDX., BDEG.) {
	# Function for B-spline basis
    dx <- (XR. - XL.)/NDX.
	knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
    P <- outer(X., knots, tpower, BDEG.)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = BDEG. + 1) / (gamma(BDEG. + 1) * dx ^ BDEG.)
    B <- (-1) ^ (BDEG. + 1) * P %*% t(D)
    res <- list(B = B, knots = knots)
	res 
}
