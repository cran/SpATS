construct.fixed.part <-
function(formula, data, genotype, grandom) {
	env <- environment(formula) 
	if(inherits(formula, "character"))          
		formula <- as.formula(formula)
	
	# If genotype fixed and also included in formula (in interaction), 
	# add it to the formula (main effect) to obtain the correct design matrix
	# It is removed later on
	if(!grandom & (genotype %in% all.vars(formula))) {
		 formula <- update(formula, as.formula(paste("~", genotype, "+ .")))
	}

	mf <- model.frame(formula, data, drop.unused.levels = TRUE)
	mt <- terms(mf)   
	X <- model.matrix(mt, mf)
	
	dim <- table(attr(X,"assign"))[-1]
	names(dim) <- attr(mt, "term.labels")
	
	attr(mt, "contrast") <- attr(X,"contrast")
	attr(mt, "xlev") <- .getXlevels(mt, mf)

	# Remove "main" effect for genotype
	if(!grandom & (genotype %in% all.vars(formula))) {
		X <- X[,-(1:(dim[1]+1)),drop = FALSE]
		dim <- dim[-1]
	} else {
		X <- X[,-1, drop = FALSE]
	}
	
	# For prediction	
	# mfp <- model.frame(mt, newdata, xlev = attr(mt, "xlev"))
	# Xp <- model.matrix(mt, data = mfp, contrasts.arg = attr(mt, "contrast"))
	
	attr(dim, "random") <- attr(dim, "sparse") <- attr(dim, "spatial") <- rep(FALSE, length(dim)) 	
	#res <- list(X = X[,-1, drop = FALSE], dim = dim, terms = mt)
	res <- list(X = X, dim = dim, terms = mt)
	res	
}
