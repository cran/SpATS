print.summary.SpATS <-
function(x, which = c("dimensions", "variances", "all"), ...) {
	print.SpATS(x)
	which <- match.arg(which)
	if(which == "variances" | which == "all") {
		cat("\nVariance components:\n")
		print(x$p.table.vc, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}
	if(which == "dimensions" | which == "all") {
		cat("\nDimensions:\n")
		print(x$p.table.dim, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)	
		cat('\nType codes: F \'Fixed\'    R \'Random\'    S \'Smooth/Semiparametric\'\n\n')
	}	
	invisible(x)	
}
