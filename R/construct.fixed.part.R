construct.fixed.part <-
function(formula, data) {
	env <- environment(formula) 
	if(inherits(formula, "character"))          
		formula <- as.formula(formula)
	
	mf <- model.frame(formula, data, drop.unused.levels = TRUE)
	mt <- terms(mf)   
	X <- model.matrix(mt, mf)
	
	dim <- table(attr(X,"assign"))[-1]
	names(dim) <- attr(mt, "term.labels")
	
	attr(mt, "contrast") <- attr(X,"contrast")
	attr(mt, "xlev") <- .getXlevels(mt, mf)
	
	# For prediction	
	# mfp <- model.frame(mt, newdata, xlev = attr(mt, "xlev"))
	# Xp <- model.matrix(mt, data = mfp, contrasts.arg = attr(mt, "contrast"))
	
	attr(dim, "random") <- attr(dim, "sparse") <- attr(dim, "spatial") <- rep(FALSE, length(dim)) 	
	res <- list(X = X[,-1, drop = FALSE], dim = dim, terms = mt)
	res	
}
