obtain.nominal.dimension <- function(MM, dim, random., spatial., weights) {
	MM <- as.matrix(MM)
	MM <- MM[weights != 0,]
	dim.nom <- dim
	random <- which(random. == TRUE & spatial. == FALSE)
	fixed.pos <- create.position.indicator(dim, !random.)

	np.e <- cumsum(dim)
	np.s <- np.e - dim + 1

	X <- MM[,fixed.pos]

	for(i in random) {
		Z <- MM[,np.s[i]:np.e[i]]
		#tst <- cbind(rowSums(Z), X)
		dim.nom[i] =  qr(cbind(Z,X))$rank - qr(X)$rank #dim[i] - (ncol(tst) - qr(tst)$rank)
	}
	dim.nom
}