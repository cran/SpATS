obtain.spatialtrend <-
function(object, grid = c(100,100), ...) {
	terms.formula <- object$terms$spatial$terms.formula
	
	grid = if(length(grid) == 1) rep(grid, 2) else grid
	
	x.coord <- object$data[,terms.formula$x.coord]
	y.coord <- object$data[,terms.formula$y.coord]

	col.p <- seq(min(x.coord), max(x.coord), l = grid[1])
	row.p <- seq(min(y.coord), max(y.coord), l = grid[2])

	B1p <- spline.bbase(object$terms$spatial$MM$MM1$knots, col.p, terms.formula$degree[1])
	B2p <- spline.bbase(object$terms$spatial$MM$MM2$knots, row.p, terms.formula$degree[2])

	X1p <- B1p%*%object$terms$spatial$MM$MM1$U.X
	X2p <- B2p%*%object$terms$spatial$MM$MM2$U.X

	Z1p <- B1p%*%object$terms$spatial$MM$MM1$U.Z
	Z2p <- B2p%*%object$terms$spatial$MM$MM2$U.Z
	
	# Coefficients associated to the spatial component
	fixed.spat.coef <- object$coeff[object$terms$spatial$fixed$pos]
	random.spat.coef <- object$coeff[object$terms$spatial$random$pos]
	
	if(terms.formula$type == "SAP") {
		Xp <- X2p%x%X1p
		Xp <- Xp[,-1,drop = FALSE]
		Zp = cbind(X2p%x%Z1p, Z2p%x%X1p, Z2p%x%Z1p)
	} else {
		B1pn <- spline.bbase(object$terms$spatial$MMn$MM1$knots, col.p, terms.formula$degree[1])
		B2pn <- spline.bbase(object$terms$spatial$MMn$MM2$knots, row.p, terms.formula$degree[2])
		Z1pn <- B1pn%*%object$terms$spatial$MMn$MM1$U.Z
		Z2pn <- B2pn%*%object$terms$spatial$MMn$MM2$U.Z

		Xp <- X2p%x%X1p
		Xp <- Xp[,-1,drop = FALSE]		
		
		# Separate for each PS-ANOVA component
		Zp1 <- X2p[,1, drop = FALSE]%x%Z1p
		Zp2 <- Z2p%x%X1p[,1, drop = FALSE]
		Zp3 <- X2p[,-1, drop = FALSE]%x%Z1p
		Zp4 <- Z2p%x%X1p[,-1, drop = FALSE]
		Zp5 <- Z2pn%x%Z1pn
		
		Zp = cbind(Zp1, Zp2, Zp3, Zp4, Zp5)
		dims <- c(ncol(Zp1), ncol(Zp2), ncol(Zp3), ncol(Zp4), ncol(Zp5))
		e <- cumsum(dims)
		s <- e - dims + 1
		# Main effectts
		eta1 <- matrix(Zp1%*%random.spat.coef[s[1]:e[1]], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		eta2 <- matrix(Zp2%*%random.spat.coef[s[2]:e[2]], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		
		# Linear-by-smooth
		eta3 <- matrix(Zp3%*%random.spat.coef[s[3]:e[3]], ncol = length(col.p), byrow = TRUE)		
		eta4 <- matrix(Zp4%*%random.spat.coef[s[4]:e[4]], ncol = length(col.p), byrow = TRUE)
		
		# Smooth-by-smooth
		eta5 <- matrix(Zp5%*%random.spat.coef[s[5]:e[5]], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
	
	}
	eta <- cbind(Xp,Zp)%*%c(fixed.spat.coef, random.spat.coef)
	res <- list(col.p = col.p, row.p = row.p, fit = matrix(eta, nrow = length(row.p), ncol = length(col.p), byrow = TRUE))
	if (terms.formula$type != "SAP")
		res$pfit <- list(fv = eta1, fu = eta2, uhv = eta3, vhu = eta4, fuv = eta5)
	res
}
