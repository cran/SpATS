construct.genotype.random.part <-
function(genotype, geno.decomp = NULL, data, weights) {
	geno <- data[, genotype]
	geno = as.factor(geno)
	decomp <- if(!is.null(geno.decomp)) data[, geno.decomp] else NULL
	
	if(!is.null(decomp)) {
		ord <- order(decomp)
		geno_names <- as.character(unique(na.omit(geno[ord])))
			
		# Genotypes presented in the sample
		geno.groups.dim <- unlist(lapply(split(data.frame(geno, weights), decomp), function(x) {
				temp <- droplevels(x[x$weights != 0, 1])
				res <- nlevels(temp)
				res}))
		names(geno.groups.dim) <- paste(geno.decomp, levels(decomp), sep = "")
		geno.groups.dim <- geno.groups.dim[geno.groups.dim != 0]
		
		# All genotypes
		pop.groups.dim <- unlist(lapply(split(data.frame(geno, weights), decomp), function(x) {
				temp <- droplevels(x[,1])
				res <- nlevels(temp)
				res}))						
		pop_names <- rep(levels(decomp)[pop.groups.dim != 0], pop.groups.dim[pop.groups.dim != 0])
		pop.groups.dim <- pop.groups.dim[pop.groups.dim != 0]	
			
	} else {
		geno_names = levels(geno)
		geno.groups.dim <- nlevels(droplevels(geno[weights != 0]))
		names(geno.groups.dim) = genotype
		pop_names <- rep("POP0", length(geno_names))
	}
	ndx = get_ndx_scores(weights = weights, geno = geno, names = geno_names)
	Z_geno = make_design_matrix(geno, geno_names)	

	# Only those genotypes presented in the dataset
	Z_geno <- Z_geno[,ndx]
	attr(Z_geno, "colnames") <- geno_names[ndx]
	geno_dim <- rep(0, length(geno_names))
	geno_dim[ndx] <- colSums(Z_geno)

	e <- cumsum(geno.groups.dim)
	s <- e - geno.groups.dim + 1

	g <- list()
	for(i in 1:length(geno.groups.dim)) {
		g[[i]] <- rep(0, sum(geno.groups.dim))
		g[[i]][s[i]:e[i]] <- 1
	}
	names(g) <- names(geno.groups.dim)
	dim <- geno.groups.dim
	attr(dim, "random") <- attr(dim, "sparse") <- rep(TRUE, length(dim))
	attr(dim, "spatial") <- rep(FALSE, length(dim))
	
	# Initialize variance components
	init.var <- rep(1, length(g))
	
	res <- list(Z = Z_geno, dim = dim, g = g, init.var = init.var, terms = list(geno_names = geno_names, geno_dim = geno_dim, pop_names = pop_names, ndx = ndx))
}
