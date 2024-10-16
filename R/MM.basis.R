MM.basis <-
function (x, xl, xr, ndx, bdeg, pord, decom = 1) {
	# Check that ndx, bdf and pord are integers
	ndx.new <- round(ndx)
	if (all.equal(ndx.new, ndx) != TRUE) {
		warning("argument 'nseg' of PSANOVA()/SAP() should be integer and has been rounded")
	}
	ndx <- ndx.new
	
	pord.new <-round(pord)
	if (all.equal(pord.new, pord) != TRUE) {
		warning("argument 'pord' of PSANOVA()/SAP() should be integer and has been rounded")
	}
	pord <- pord.new

	bdeg.new <-round(bdeg)
	if (all.equal(bdeg.new, bdeg) != TRUE) {
		warning("argument 'degree' of PSANOVA()/SAP() should be integer and has been rounded")
	}
	bdeg <- bdeg.new


	Bb = bbase(x,xl,xr,ndx,bdeg)
	knots <- Bb$knots
	B = Bb$B
	m = ncol(B)
	n = nrow(B)
	D = diff(diag(m), differences=pord)
	P.svd = svd(crossprod(D))
	U.Z = (P.svd$u)[,1:(m-pord)] # eigenvectors
	d = (P.svd$d)[1:(m-pord)]  # eigenvalues
	Z = B%*%U.Z
	U.X = NULL
	if(decom == 1) {
		U.X = ((P.svd$u)[,-(1:(m-pord))])
		X = B%*%U.X
	} else if (decom == 2){
		X = NULL
		for(i in 0:(pord-1)){
			X = cbind(X,x^i)
		}
	} else if(decom == 3) {
		U.X = NULL
		for(i in 0:(pord-1)){
			U.X = cbind(U.X,knots[-c((1:pord),(length(knots)- pord + 1):length(knots))]^i)
		}
		X = B%*%U.X
	} else if(decom == 4) { # Wood's 2013
		X = B%*%((P.svd$u)[,-(1:(m-pord))])
		id.v <- rep(1, nrow(X))
		D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
		Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
		X <- X%*%Xf
		U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
	}
	list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
}
