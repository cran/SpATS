construct.fixed.prediction.matrix <-
function(object, newdata) {
	if(!is.null(object$terms$fixed)) {
		mfp <- model.frame(object$terms$fixed, newdata, xlev = attr(object$terms$fixed, "xlev"))
		Xp <- model.matrix(object$terms$fixed, data = mfp, contrasts.arg = attr(object$terms$fixed, "contrast"))
		Xp <- Xp[,-1,drop = FALSE]
	} else {
		Xp <- NULL
	}
	Xp
}
