predict.SpATS <-
function(object, newdata = NULL, which = NULL, ...) {
	if((is.null(newdata) & is.null(which)) | (!is.null(newdata) & !is.null(which))) {
		stop("Either newdata or which must be indicated")
	}
	if(!is.null(which) & !is.character(which))
		stop("Which argument must be a character vector")
		
	if(!is.null(newdata) & !is.data.frame(newdata))
		stop("newdata argument must be a data frame")
			
	model.terms <- c(object$model$geno$genotype, 
					 object$terms$spatial$terms.formula$x.coord, 
					 object$terms$spatial$terms.formula$y.coord,
					 if(!is.null(object$terms$fixed)) all.vars(object$terms$fixed),
					 if(!is.null(object$terms$random)) all.vars(object$terms$random))
	if(!is.null(newdata)) {
		if(any(!(model.terms %in% names(newdata))))
			stop("Not all needed variables are supplied in newdata")
	}
	if(!is.null(which)) {
		which <- which[which %in% model.terms]
		if((object$terms$spatial$terms.formula$x.coord %in% which & !(object$terms$spatial$terms.formula$y.coord%in%which)) | (object$terms$spatial$terms.formula$y.coord %in% which & !(object$terms$spatial$terms.formula$x.coord %in% which)))
			stop(paste("Both ",object$terms$spatial$terms.formula$x.coord," and ", object$terms$spatial$terms.formula$y.coord," must be supplied in argument which", sep = ""))
		
		newdata <- object$data[,which,drop=FALSE]
		newdata <- newdata[!duplicated(newdata),,drop = FALSE]
		
		# Get only the genotypes in the sample
		if(object$model$geno$genotype %in% which) {
			geno.in.sample <- object$terms$geno$geno_names[object$terms$geno$ndx]
			ind.geno <- newdata[,object$model$geno$genotype] %in% geno.in.sample
			newdata <- newdata[ind.geno,,drop = FALSE]
		}

		for(i in model.terms[!model.terms %in% which]) {
			if(is.factor(object$data[,i])) {
				if(!is.null(object$terms$fixed) & (i %in% all.vars(object$terms$fixed))) {
					newdata[,i] <- attr(object$terms$fixed, "xlev")[[i]][1]
				} else if (i %in% object$model$geno$genotype & !object$model$geno$as.random) {
					newdata[,i] <- as.factor(object$terms$geno$geno_names[object$terms$geno$ndx[1]]) # Reference
				} else {
					newdata[,i] <- as.factor(NA)
				}
			} else {
				newdata[,i] <- mean(object$data[,i], na.rm = TRUE)
			}
		}
		newdata <- newdata[order(newdata[,which]),]						
	}
	# NAs in the covariates: allowed in the random part, as well as the genotype when random
	na.terms <- c(object$terms$spatial$terms.formula$x.coord, 
					 object$terms$spatial$terms.formula$y.coord,
					 if(!object$model$geno$as.random) object$model$geno$genotype,
					 if(!is.null(object$terms$fixed)) all.vars(object$terms$fixed))
					 
	na.ind <- apply(is.na(newdata[,na.terms]), 1, any)
	newdata <- newdata[!na.ind,]
	
	# Construct matrices
	spatial.part <- construct.spatial.prediction.matrix(object, newdata)
	if(!is.null(which)) {
		# Delete spatial prediction if the spatial coordinates have not been specified in the argument 'which'
		if(!(object$terms$spatial$terms.formula$x.coord %in% which) & !(object$terms$spatial$terms.formula$y.coord%in%which)) {
			spatial.part$X <- matrix(0, ncol = ncol(spatial.part$X), nrow = nrow(spatial.part$X))
			spatial.part$Z <- matrix(0, ncol = ncol(spatial.part$Z), nrow = nrow(spatial.part$Z))
			newdata[,object$terms$spatial$terms.formula$x.coord] <- "Excluded"
			newdata[,object$terms$spatial$terms.formula$y.coord] <- "Excluded"
		}
	}
	MMsp <- construct.genotype.prediction.matrix(object, newdata)
	MMnsp <- cbind("Intercept" = rep(1, nrow(newdata)), construct.fixed.prediction.matrix(object, newdata), construct.random.prediction.matrix(object, newdata), spatial.part$X, spatial.part$Z)
	
	# Predicted values
	predicted.values <- as.vector(cbind.spam(MMsp, MMnsp)%*%as.vector(object$coeff))
	# Standard errors
	standard.errors <- sqrt(rowSums(MMsp%*%object$vcov$C11_inv*MMsp) + rowSums(2*MMsp%*%object$vcov$C12_inv*MMnsp) + rowSums(MMnsp%*%object$vcov$C22_inv*MMnsp))
	
	na.terms <- c(if(object$model$geno$as.random) object$model$geno$genotype, if(!is.null(object$terms$random)) all.vars(object$terms$random))
	if(!is.null(na.terms)) {
		for(i in na.terms) {
			if(sum(is.na(newdata[,i])) == nrow(newdata)) {
				levels(newdata[,i]) <- c(levels(newdata[,i]), "Excluded")
				newdata[,i] <- "Excluded"
			}
				
		}
	}
	res <- newdata
	res$predicted.values <- predicted.values
	res$standard.errors <- standard.errors
	res	
}