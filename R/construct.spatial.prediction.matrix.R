construct.spatial.prediction.matrix <-
function(object, newdata) {
	terms.formula <- object$terms$spatial$terms.formula
	
	x.coord <- newdata[,terms.formula$x.coord]
	y.coord <- newdata[,terms.formula$y.coord]

	B1p <- spline.bbase(object$terms$spatial$MM$MM1$knots, x.coord, terms.formula$degree[1])
	B2p <- spline.bbase(object$terms$spatial$MM$MM2$knots, y.coord, terms.formula$degree[2])

	X1p <- B1p%*%object$terms$spatial$MM$MM1$U.X
	X2p <- B2p%*%object$terms$spatial$MM$MM2$U.X

	Z1p <- B1p%*%object$terms$spatial$MM$MM1$U.Z
	Z2p <- B2p%*%object$terms$spatial$MM$MM2$U.Z

	Xp = Rten2(X2p, X1p)
	Xp <- Xp[,-1,drop = FALSE]
	
	# Nested bases
	B1pn <- spline.bbase(object$terms$spatial$MMn$MM1$knots, x.coord, terms.formula$degree[1])
	B2pn <- spline.bbase(object$terms$spatial$MMn$MM2$knots, y.coord, terms.formula$degree[2])
	
	Z1pn <- B1pn%*%object$terms$spatial$MMn$MM1$U.Z
	Z2pn <- B2pn%*%object$terms$spatial$MMn$MM2$U.Z
				
	if(terms.formula$type == "SAP") {
		Zp = cbind(Rten2(X2p, Z1p), Rten2(Z2p, X1p), Rten2(Z2pn, Z1pn))
	} else {
		Zp <- cbind(Rten2(X2p[,1, drop = FALSE], Z1p), Rten2(Z2p, X1p[,1, drop = FALSE]), Rten2(X2p[,-1, drop = FALSE], Z1p), Rten2(Z2p, X1p[,-1, drop = FALSE]), Rten2(Z2pn, Z1pn))
	}
	res <- list(X = Xp, Z = Zp)
	res	
}
