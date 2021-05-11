construct.design.matrix <-
function(genotype, geno.decomp = NULL, grandom = FALSE, spatial, fixed = NULL, random = NULL, data, weights, na.res) {
	# Create the mixed model matrices
	MMs <- MMns <- NULL
	gg <- init.var <- dim <- list()
		
	# Genotype
	if(grandom) {
		geno.part <- construct.genotype.random.part(genotype = genotype, geno.decomp = geno.decomp, data = data, weights = weights)
		MMs <- cbind.spam(MMs, geno.part$Z)
		dim <- c(dim, list(geno.part$dim))
		gg <- c(gg, list(geno.part$g))
		init.var <- c(init.var, list(geno.part$init.var))
	} else {
		geno.part <- construct.genotype.fixed.part(genotype = genotype, data = data, weights = weights)
		MMs <- cbind.spam(MMs, geno.part$X)
		dim <- c(dim, list(geno.part$dim))
	}
	# Fixed part
	# Intercept
	int <- rep(1, nrow(data))
	dim.int <- c("Intercept" = 1)
	attr(dim.int, "random") <- FALSE
	attr(dim.int, "spatial") <- FALSE
	attr(dim.int, "sparse") <- FALSE	
	MMns <- cbind(MMns, "Intercept" = int)
	dim <- c(dim, list(dim.int))
	
	if(!is.null(fixed)) {
		fixed.part <- construct.fixed.part(formula = fixed, data = data)
		MMns <- cbind(MMns, fixed.part$X)
		dim <- c(dim, list(fixed.part$dim))
	}
	# Random part
	if(!is.null(random)) {
		random.part <- construct.random.part(formula = random, data = data)
		MMns <- cbind(MMns, random.part$Z)
		dim <- c(dim, list(random.part$dim))
		gg <- c(gg, list(random.part$g))
		init.var <- c(init.var, list(random.part$init.var))
	}
	# Smooth (spatial part)
	spat.part <- construct.2d.pspline(formula = spatial, data = data, na.res = na.res)
	MMns <- cbind(MMns, spat.part$X, spat.part$Z)
	dim <- c(dim, list(spat.part$dim$fixed), list(spat.part$dim$random))
	gg <- c(gg, list(spat.part$g))
	init.var <- c(init.var, list(spat.part$init.var))

	# Capital lambda
	g <- construct.capital.lambda(gg)

	# Type effect indicators
	random. <- unlist(lapply(dim, get.attribute, "random"))
	spatial. <- unlist(lapply(dim, get.attribute, "spatial"))

	spat.part$terms$fixed$pos <- create.position.indicator(unlist(dim), !random. & spatial.)
	spat.part$terms$random$pos <- create.position.indicator(unlist(dim), random. & spatial.)
	
	res <- list()
	res$geno.part <- geno.part
	res$fixed.part <- if(!is.null(fixed)) fixed.part else NULL
	res$random.part <- if(!is.null(random)) random.part else NULL
	res$spat.part <- spat.part	
	# Terms	
	res$terms$spatial <- spat.part$terms
	res$terms$fixed <- if(!is.null(fixed)) fixed.part$terms else NULL
	res$terms$random <- if(!is.null(random)) random.part$terms else NULL
	res$terms$geno <- geno.part$terms
	# Capital Lambda	
	res$g <- g
	# Dimension
	res$dim <- dim
	# Intial values for the variance compoments
	res$init.var <- init.var
	# Design Matrices
	res$MM <- list(MMs = MMs, MMns = MMns)
	res
}
