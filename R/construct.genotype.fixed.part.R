construct.genotype.fixed.part <-
function(genotype, data, weights) {
	geno <- data[,genotype]
	geno <- as.factor(geno)
	geno_names = levels(geno)
	X_geno_c = make_design_matrix(geno, geno_names)  
	ndx = get_ndx_scores(weights, geno, geno_names)
	geno_dim <- colSums(X_geno_c)
	
	# Only those genotypes presented in the dataset (minus the first for identifiability reasons)
	geno_dim <- rep(0, length(geno_names))
	geno_dim[ndx] <- colSums(X_geno_c[,ndx])
	X_geno = X_geno_c[, ndx[2:length(ndx)], drop = FALSE]
	
	attr(X_geno, "colnames") <- geno_names[ndx[2:length(ndx)]]
	
	dim <- ncol(X_geno)
	names(dim) <- genotype
	
	attr(dim, "random") <- attr(dim, "spatial") <- rep(FALSE, length(dim)) 
	attr(dim, "sparse") <- rep(TRUE, length(dim))
	
	res <- list(X = X_geno, dim = dim, terms  = list(geno_names = geno_names, geno_dim = geno_dim, ndx = ndx))
}
