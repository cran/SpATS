print.SpATS <-
function(x, ...) {	
	genotype.lab <- ifelse(x$model$geno$as.random, "Genotypes (as random):", "Genotypes (as fixed):")
	text <- paste('\nSpatial analysis of trials with splines \n\n', 
				  sprintf('%-27s %-10s', 'Response:', x$model$response), '\n',
				  sprintf('%-27s %-10s', genotype.lab, x$model$geno$genotype), '\n',
				  sprintf('%-27s %-10s', 'Spatial:', Reduce(paste, deparse(x$model$spatial, width.cutoff = 200L))),
				  if(!is.null(x$model$fixed)) paste('\n', sprintf('%-27s %-10s', 'Fixed:', Reduce(paste, deparse(x$model$fixed, width.cutoff = 200L))), sep = '') else NULL,
				  if(!is.null(x$model$random)) paste('\n', sprintf('%-27s %-10s', 'Random:', Reduce(paste, deparse(x$model$random, width.cutoff = 200L))), sep = '') else NULL,
				  '\n\n\n',
				  sprintf('%-30s %.0f', 'Number of observations:', x$nobs),'\n',
				  sprintf('%-30s %.0f', 'Number of missing data:', nrow(x$data) - x$nobs), '\n',
				  sprintf('%-30s %.2f', 'Effective dimension:', sum(x$eff.dim)), '\n',
				  sprintf('%-30s %.3f', 'Deviance:', x$dev), '\n', sep = '')
	cat(text)
	invisible(x)
}
