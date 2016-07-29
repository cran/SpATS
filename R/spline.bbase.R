spline.bbase <-
function(knots, X., BDEG.) {
	dx <- diff(knots)[1]
	P <- outer(X., knots, tpower, BDEG.)
	n <- dim(P)[2]
	D <- diff(diag(n), diff = BDEG. + 1)/(gamma(BDEG. + 1)*dx^BDEG.)
	B <- (-1) ^ (BDEG. + 1) * P %*% t(D)
	B
}
