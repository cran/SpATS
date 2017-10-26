construct.henderson.matrix <-
function(mat, la, Ginv, dim, sparse., random., fixed.matrices = NULL, as.random = FALSE) {
	fixed.pos.ns <- create.position.indicator(dim[!sparse.], !random.[!sparse.])
	sparse.random.pos <- create.position.indicator(dim, random. & sparse.)	
	Ginvns <- rep(0, l = sum(dim[!sparse.]))	
	if(as.random) {
		Ginvs <- Ginv[sparse.random.pos]
		Ginvns[-fixed.pos.ns] <- Ginv[-sparse.random.pos]		
		C11 <- (1/la[1])*mat$XtX + spam::diag.spam(Ginvs)
		C12 <- (1/la[1])*mat$XtZ.
		C21 <- (1/la[1])*mat$ZtX.
		C22 <- (1/la[1])*mat$ZtZ. + diag(Ginvns)
	} else {
		Ginvns[-fixed.pos.ns] <- Ginv
		C11 <- (1/la[1])*mat$XtX.
		C12 <- (1/la[1])*mat$XtZ.
		C22 <- (1/la[1])*t(mat$ZtZ.) + diag(Ginvns)
		C21 <- (1/la[1])*mat$ZtX.
	}
	A = C11
	A_inv = 1.0/A
	if(is.null(fixed.matrices) | is.null(fixed.matrices$M))
		M = t(A_inv%*%t(C21))
	else
		M = fixed.matrices$M
		
	if(is.null(fixed.matrices) | is.null(fixed.matrices$C21_C11_inv_C12))
		K <- C22 - M%*%C12
	else
		K <- C22 - (1/la[1])*fixed.matrices$C21_C11_inv_C12
	
	if(is.null(fixed.matrices) | is.null(fixed.matrices$M_Xty.))
		M_Xty. <- M%*%mat$Xty.
	else
		M_Xty. <- fixed.matrices$M_Xty.
	
	res <- list()
	res$C11 <- C11
	res$C12 <- C12
	res$C21 <- C21
	res$C22 <- C22
	res$A_inv <- A_inv
	res$M <- M
	res$K <- K
	res$M_Xty. <- M_Xty.
	res
}
