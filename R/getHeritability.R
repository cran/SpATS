getHeritability <- function(object, ...) {
	if(!object$model$geno$as.random)
		stop("Heritability can only be calculated when genotype is random")
	if(!is.null(object$model$geno$geno.decomp)) {
		geno.decomp <- object$model$geno$geno.decomp
		decomp <- unique(object$terms$geno$pop_names)
		select <- paste(geno.decomp, decomp, sep = "")
		dim <- object$dim[select]
	} else {
		select <- object$model$geno$genotype
		dim <- object$dim[select]
	}

	ed.geno <- object$eff.dim[select]/(dim - 1)
	names(ed.geno) <- select
	res <- round(ed.geno, 2)
	res
}