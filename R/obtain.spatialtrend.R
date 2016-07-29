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
		Zp = cbind(X2p[,1]%x%Z1p, Z2p%x%X1p[,1], X2p[,-1]%x%Z1p, Z2p%x%X1p[,-1], Z2pn%x%Z1pn)
		Xp <- Xp[,-1,drop = FALSE]
	}
	p.fit <- cbind(Xp, Zp)%*%object$coeff[c(object$terms$spatial$fixed$pos, object$terms$spatial$random$pos)]
	res <- list(col.p = col.p, row.p = row.p, fit = matrix(p.fit, nrow = length(row.p), ncol = length(col.p), byrow = TRUE))
	res
}
