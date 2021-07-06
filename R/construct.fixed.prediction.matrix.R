construct.fixed.prediction.matrix <-
function(object, newdata, genotype, grandom) {
	if(!is.null(object$terms$fixed)) {
		mfp <- model.frame(object$terms$fixed, newdata, xlev = attr(object$terms$fixed, "xlev"))
		Xp <- model.matrix(object$terms$fixed, data = mfp, contrasts.arg = attr(object$terms$fixed, "contrast"))

		if(!grandom & (genotype %in% all.vars(object$terms$fixed))) {
			dim <- table(attr(Xp,"assign"))[-1]
			Xp <- Xp[,-(1:(dim[1]+1)),drop = FALSE]
		} else {
			Xp <- Xp[,-1,drop = FALSE]
		}
	} else {
		Xp <- NULL
	}
	Xp
}
