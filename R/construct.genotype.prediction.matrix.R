construct.genotype.prediction.matrix <-
function(object, newdata) {
	Z_geno <- make_design_matrix(newdata[,object$model$geno$genotype], object$terms$geno$geno_names)
	if(object$model$geno$as.random)
		Z_geno <- Z_geno[,object$terms$geno$ndx, drop = FALSE]
	else
		Z_geno <- Z_geno[, object$terms$geno$ndx[2:length(object$terms$geno$ndx)], drop = FALSE]
	Z_geno	
}
