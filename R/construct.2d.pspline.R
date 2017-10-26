construct.2d.pspline <-
function(formula, data) {
	env <- environment(formula) 
	if(inherits(formula, "character"))          
		formula <- as.formula(formula)
		
	res <- interpret.SpATS.formula(formula)
	
	x1 <- data[ ,res$x.coord]
	x2 <- data[ ,res$y.coord]
	
	type = res$type
	
	MM1 = MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 4)
	MM2 = MM.basis(x2, min(x2), max(x2), res$nseg[2], res$degree[2], res$pord[2], 4)
		
	X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
	X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B
	
	c1 = ncol(B1); c2 = ncol(B2)
	
	# Nested bases
	if(res$nest.div[1] == 1) {
		MM1n <- MM1
		Z1n <- Z1
		c1n <- c1
		d1n <- d1	
	} else {
		MM1n = MM.basis(x1, min(x1), max(x1), res$nseg[1]/res$nest.div[1], res$degree[1], res$pord[1], 4)
		Z1n <- MM1n$Z
		d1n <- MM1n$d
		c1n <-  ncol(MM1n$B)  					
	}
	if(res$nest.div[2] == 1) {
		MM2n <- MM2
		Z2n <- Z2
		c2n <- c2
		d2n <- d2	
	} else {
		MM2n = MM.basis(x2, min(x2), max(x2), res$nseg[2]/res$nest.div[2], res$degree[2], res$pord[2], 4)
		Z2n <- MM2n$Z
		d2n <- MM2n$d
		c2n <-  ncol(MM2n$B)  					
	}
	
	x.fixed <- y.fixed <- ""
	for(i in 0:(res$pord[1]-1)){
		if(i == 1) 
			x.fixed <- c(x.fixed, res$x.coord)
		else if( i > 1)
			x.fixed <- c(x.fixed, paste(res$x.coord, "^", i, sep = ""))
	}
	for(i in 0:(res$pord[2]-1)){
		if(i == 1) 
			y.fixed <- c(y.fixed, res$y.coord)
		else if( i > 1)
			y.fixed <- c(y.fixed, paste(res$y.coord, "^", i, sep = ""))
	}
	xy.fixed <- NULL
	for(i in 1:length(y.fixed)) {
		xy.fixed <- c(xy.fixed, paste(y.fixed[i], x.fixed, sep= ""))
	}
	xy.fixed <- xy.fixed[xy.fixed != ""]
	names.fixed <- xy.fixed
	
	smooth.comp <- paste("f(", res$x.coord,",", res$y.coord,")", sep = "")
	
	if(type == "SAP") {
		names.random <- paste(smooth.comp, c(res$x.coord, res$y.coord), sep = "|")				
		X = Rten2(X2, X1)		
		# Delete the intercept
		X <- X[,-1,drop = FALSE]
		Z = cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2n, Z1n))
		
		dim.random <- c((c1 -res$pord[1])*res$pord[2] , (c2 - res$pord[2])*res$pord[1], (c1n - res$pord[1])*(c2n - res$pord[2]))		
		dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
		names(dim$fixed) <- names.fixed
		names(dim$random) <- paste(smooth.comp, "Global")
		
		# Variance/Covariance components
		g1u <- rep(1, res$pord[2])%x%d1
		g2u <- d2%x%rep(1, res$pord[1])
		g1b <- rep(1, c2n - res$pord[2])%x%d1n
		g2b <- d2n%x%rep(1, c1n - res$pord[1])
		
		g <- list()	
		g[[1]] <- c(g1u, rep(0, dim.random[2]), g1b)
		g[[2]] <- c(rep(0, dim.random[1]), g2u, g2b)
		
		names(g) <- names.random
		
	} else {		
		one1. <- X1[,1, drop = FALSE]
		one2. <- X2[,1, drop = FALSE]
		
		x1. <- X1[,-1, drop = FALSE]
		x2. <- X2[,-1, drop = FALSE]
		
		# Fixed and random matrices
		X <- Rten2(X2, X1)
		# Delete the intercept
		X <- X[,-1,drop = FALSE]
		Z <- cbind(Rten2(one2., Z1), Rten2(Z2, one1.), Rten2(x2., Z1), Rten2(Z2, x1.), Rten2(Z2n, Z1n))
		
		dim.random <- c((c1-res$pord[1]), (c2-res$pord[2]), (c1-res$pord[1])*(res$pord[2]-1), (c2-res$pord[2])*(res$pord[1]-1), (c1n-res$pord[2])*(c2n-res$pord[2]))
			
		# Variance/Covariance components		
		g1u <- d1
		g2u <- d2
		
		g1v <- rep(1, res$pord[2] - 1)%x%d1
		g2v <- d2%x%rep(1,res$pord[1] - 1)
		
		g1b <- rep(1, c2n - res$pord[2])%x%d1n
		g2b <- d2n%x%rep(1, c1n - res$pord[1])
		
		g <- list()
		
		if(type == "SAP.ANOVA") {
			g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
			g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
			g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, dim.random[4]), g1b)
			g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, g2b)
	
			names.random <- c(paste("f(", res$x.coord,")", sep = ""), paste("f(", res$y.coord,")", sep = ""), paste(smooth.comp, c(res$x.coord, res$y.coord), sep = "|"))			
			dim <- list(fixed = rep(1, ncol(X)), random = c(dim.random[1:2], sum(dim.random[-(1:2)])))		
			names(dim$fixed) <- names.fixed
			names(dim$random) <- c(names.random[1:2], paste(smooth.comp, "Global"))
			names(g) <- names.random
		} else {
			g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
			g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
			g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, sum(dim.random[4:5])))
			g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, rep(0, dim.random[5]))
			g[[5]] <- c(rep(0, sum(dim.random[1:4])), g1b + g2b)
			
			names.random <- c(paste("f(", res$x.coord,")", sep = ""), paste("f(", res$y.coord,")", sep = ""),
							paste("f(", res$x.coord,"):", res$y.coord, sep = ""),
							paste(res$x.coord,":f(", res$y.coord,")", sep = ""),
							paste("f(", res$x.coord,"):f(", res$y.coord,")", sep = ""))
							
			dim <- list(fixed = rep(1, ncol(X)), random = dim.random)		
			names(dim$fixed) <- names.fixed
			names(dim$random) <- names.random
			names(g) <- names.random
		}		
	}
	colnames(X) <- names.fixed
	colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
	
	attr(dim$fixed, "random") <- attr(dim$fixed, "sparse") <- rep(FALSE, length(dim$fixed))
	attr(dim$fixed, "spatial") <- rep(TRUE, length(dim$fixed))
	
	attr(dim$random, "random") <- attr(dim$random, "spatial") <- rep(TRUE, length(dim$random)) 
	attr(dim$random, "sparse") <- rep(FALSE, length(dim$random))

	terms <- list()
	terms$MM <- list(MM1 = MM1, MM2 = MM2)
	terms$MMn <- list(MM1 = MM1n, MM2 = MM2n)
	terms$terms.formula <- res
	
	attr(terms, "term") <- smooth.comp
	
	# Initialize variance components
	init.var <- rep(1, length(g))
	
	res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var, terms = terms)	
	res
}
