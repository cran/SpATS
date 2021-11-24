construct.matrices <-
function(X, Z, z, w) {
	XtW. = t(sweep(X, MARGIN = 1, w, '*')) #XtW. = t(X*w)
	XtX. = XtW.%*%X
	XtZ. = XtW.%*%Z
	ZtX. = t(XtZ.)
	ZtW. =  t(Z*w)
	ZtZ. = ZtW.%*%Z
	Xty. = XtW.%*%z
	Zty. = ZtW.%*%z
	yty. <- sum((z^2)*w)
	ZtXtZ = rbind(XtZ., ZtZ.)
	u <- c(Xty.,Zty.)
	res <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
}
