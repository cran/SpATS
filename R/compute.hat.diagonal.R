compute.hat.diagonal <-
function(mat, la, A_inv, K_inv, M, Ginv, dim, sparse., random., as.random = FALSE) {
	fixed.pos.ns <- create.position.indicator(dim[!sparse.], !random.[!sparse.])
	sparse.random.pos <- create.position.indicator(dim, random. & sparse.)
	Ginvns <- rep(1, l = sum(dim[!sparse.]))
	if(as.random) {
		Ginvs <- Ginv[sparse.random.pos]	
		Ginvns[-fixed.pos.ns] <- Ginv[-sparse.random.pos]
	} else {
		Ginvns[-fixed.pos.ns] <- Ginv
	}
	C21_inv <- -K_inv%*%M
	C12_inv <- t(C21_inv)
	C22_inv <- K_inv
	Hns <- (cbind(C21_inv*Ginvns, C22_inv*Ginvns))[-fixed.pos.ns,]
	dZtPZ.ns <- 1/la[1]*colSums((t(Hns)*mat$ZtXtZ[,-fixed.pos.ns]))
	if(as.random) {
		C11_inv_diag <- diag(A_inv) + rowSums(-C12_inv*t(M))
		dZtPZ.s <-  diag(C11_inv_diag*Ginvs*mat$XtX.) + colSums(t(C12_inv*Ginvs)*mat$ZtX.)
		dZtPZ <- c(1/la[1]*dZtPZ.s, dZtPZ.ns)
	} else {
		C11_inv_diag <- NULL
		dZtPZ = dZtPZ.ns
	}
	res <- list()
	res$dZtPZ <- dZtPZ
	res$C11_inv_diag <- C11_inv_diag
	res$C12_inv <- C12_inv
	res$C21_inv <- C21_inv
	res$C22_inv <- C22_inv
	res	
}
