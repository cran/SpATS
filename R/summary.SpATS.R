summary.SpATS <-
function(object, which = c("dimensions", "variances", "all"), ...) {
	which <- match.arg(which)
	# Variance components
	var.comp <- object$var.comp
	psi <- object$psi[1]
	nterms <- length(var.comp)
	model <- names(var.comp)
	col.names <- c("Variance", "SD", "log10(lambda)")
	row.names <- c(model, NA, "Residual")
	vc <- matrix(ncol = 3, nrow = nterms + 2, dimnames = list(row.names, col.names))
	vc[,1] <- c(sprintf("%.3e", var.comp), NA, sprintf("%.3e", psi))
	vc[,2] <- c(sprintf("%.3e", sqrt(var.comp)), NA, sprintf("%.3e", sqrt(psi)))
	vc[,3] <- c(sprintf("%.5f", log10(psi/var.comp)), NA, NA)
	
	# Dimensions	
	eff.dim <- object$eff.dim
	dim <- object$dim

	tot_ed <- sum(eff.dim, na.rm = TRUE)
	tot_dim <- sum(dim, na.rm = TRUE)

	dim.new <- dim[match(names(eff.dim), names(dim))]
	type <- rep(NA, length = length(dim.new))
	type[(attr(dim, "random") & !attr(dim, "spatial"))[match(names(eff.dim), names(dim))]] <- "R"
	type[(!attr(dim, "random") & !attr(dim, "spatial"))[match(names(eff.dim), names(dim))]] <- "F"
	type[is.na(type)] <- "S"
	
	eff.dim.new <- eff.dim
	smooth.comp <- attr(object$terms$spatial, "term")
	if(paste(smooth.comp, "Global") %in% names(dim)) {
		dim.new <- c(dim.new, dim[paste(smooth.comp, "Global")])
		eff.dim.new <- c(eff.dim.new, sum(eff.dim.new[grep(smooth.comp, names(eff.dim.new), fixed = TRUE)]))
		names(eff.dim.new)[length(eff.dim.new)] <- paste(smooth.comp, "Global")
		type <- c(type, "S")
	}
	# Order dimensions according to the type of effect: Fixed, random and smooth
	ord <- c(which(type == "F"), which(type == "R"), which(type == "S"))
	eff.dim.new <- eff.dim.new[ord]
	dim.new <- dim.new[ord]
	model <- model[ord]
	type <- type[ord]
		
	dim.nom <- dim.new
	dim.nom[type == "R"] <- dim.nom[type == "R"] - 1

	nterms <- length(eff.dim.new)
	model <- names(eff.dim.new)

	Nobs <- object$nobs
	
	col.names <- c("Effective", "Model", "Nominal", "Ratio", "Type")
	row.names <- c(model, NA, "Total", "Residual", "Nobs")	
	m <- matrix(ncol = 5, nrow = nterms + 4, dimnames = list(row.names,col.names))
	m[,1] <- c(sprintf("%.1f", eff.dim.new), NA, sprintf("%.1f", tot_ed), sprintf("%.1f", Nobs - tot_ed), sprintf("%.0f", Nobs))
	m[,2] <- c(sprintf("%.0f", dim.new), NA, sprintf("%.0f", tot_dim), NA, NA)
	m[,3] <- c(sprintf("%.0f", dim.nom), NA, sprintf("%.0f", sum(dim.nom, na.rm = TRUE)), NA, NA)
	m[,4] <- c(sprintf("%.2f", eff.dim.new/dim.nom), NA, sprintf("%.2f", tot_ed/sum(dim.nom, na.rm = TRUE)), NA, NA)
	m[,5] <- c(type, NA, NA, NA, NA)
	
	object$p.table.vc <- vc
	object$p.table.dim <- m
	class(object) <- "summary.SpATS"
	print(object, which = which)
}
