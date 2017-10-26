construct.random.part <-
function(formula, data) {
	env <- environment(formula) 
	if(inherits(formula, "character"))          
		formula <- as.formula(formula)
		
    mf <- model.frame(formula, data, drop.unused.levels = TRUE, na.action = NULL)
	mt <- terms(mf)    
	#f.terms <- attr(mt, "term.labels")[attr(mt,"dataClasses") == "factor"]
	f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
	Z <- model.matrix(mt, data = mf, contrasts.arg = lapply(mf[,f.terms, drop = FALSE], contrasts, contrasts = FALSE))
	Z[is.na(Z)] <- 0
	
	attr(mt, "contrast") <- attr(Z,"contrast")
	attr(mt, "xlev") <- .getXlevels(mt, mf)
	
	# For prediction	
	# mfp <- model.frame(mt, newdata, xlev = attr(mt, "xlev"))
	# Xp <- model.matrix(mt, data = mfp, contrasts.arg = attr(mt, "contrast"))
	
	dim <- table(attr(Z,"assign"))[-1]
	
	e <- cumsum(dim)
	s <- e - dim + 1

	g <- list()
	for(i in 1:length(dim)) {
		g[[i]] <- rep(0,sum(dim))
		g[[i]][s[i]:e[i]] <- 1
	}
	names(g) <- names(dim) <- attr(mt,"term.labels")
	attr(dim, "random") <- rep(TRUE, length(dim)) 
	attr(dim, "sparse") <- attr(dim, "spatial") <- rep(FALSE, length(dim))
	
	# Initialize variance components
	init.var <- rep(0.01, length(g))
		
	res <- list(Z = Z[,-1, drop = FALSE], dim = dim, g = g, init.var = init.var, terms = mt)
	res
}
