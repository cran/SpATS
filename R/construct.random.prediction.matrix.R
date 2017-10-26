construct.random.prediction.matrix <-
function(object, newdata) {
	if(!is.null(object$terms$random)) {
		Zp <- NULL
		#for(i in attr(object$terms$random,"term.labels")) {
		#	formula <- as.formula(paste("~", i, sep = ""))
		#	lev <- list()
		#	lev[[i]] <- attr(object$terms$random, "xlev")[[i]]
		#	mfp <- model.frame(formula, newdata, xlev = lev, na.action = na.pass)
		#	ctr <- list()
		#	ctr[[i]] <- attr(object$terms$random, "contrast")[[i]]
		#	temp <- model.matrix(formula, data = mfp, contrasts.arg = ctr)
		#	temp <- temp[,-1,drop = FALSE]			
		#	Zp <- cbind(Zp, temp)			
		#}
		mfp <- model.frame(object$terms$random, newdata, xlev = attr(object$terms$random, "xlev"), na.action = na.pass)
		Zp <- model.matrix(object$terms$random, data = mfp, contrasts.arg = attr(object$terms$random, "contrast"))
		Zp <- Zp[,-1,drop = FALSE]
		Zp[is.na(Zp)] <- 0
	} else {
		Zp <- NULL
	}
	Zp
}
