SpATS <-
function(response, genotype, geno.decomp = NULL, genotype.as.random = FALSE, spatial, fixed = NULL, random = NULL, data, family = gaussian(), offset = 0, weights = NULL, control = controlSpATS()) {
	control <- do.call("controlSpATS", control)
	
	if (control$monitoring) start = proc.time()[3]
	
	if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    if (inherits(data, what = 'data.frame')) {
		data <- as.data.frame(data)
	} else {
		stop("The object specified in argument 'data' is not a data frame")
	}
	
	weights <- as.vector(weights)
	if(is.null(weights)) weights = rep(1, nrow(data))
	if(length(offset) == 1) offset <- rep(offset, nrow(data))

	if(inherits(fixed, "character"))
		fixed <- as.formula(fixed)
	if(inherits(random, "character"))
		random <- as.formula(random)
	if(inherits(spatial, "character"))
		spatial <- as.formula(spatial)
	
	if(!all(sapply(all.vars(random), function(x, data) is.factor(data[,x]), data = data))) {
		stop("All variables indicated in argument 'random' should be factors")
	}

	sf <- interpret.SpATS.formula(spatial)
	
	# NAs in the covariates
	model.terms <- c(sf$x.coord, sf$y.coord, if(genotype.as.random) geno.decomp else genotype, if(!is.null(fixed)) all.vars(fixed))
	na.ind <- apply(is.na(data[,model.terms]), 1, any)
	na.pos <- (1:nrow(data))[!na.ind]
	weights <- weights*(!na.ind)*(!is.na(data[,response]))

	data.na <- data[!na.ind,]
	weights.na <- weights[!na.ind]
	offset.na <-  offset[!na.ind]
	
	y <- data.na[,response]
	nobs <- length(y[weights.na != 0])

	MM <- construct.design.matrix(genotype = genotype, geno.decomp = geno.decomp, grandom = genotype.as.random, spatial = spatial, fixed = fixed, random = random, data = data.na, weights = weights.na)
	MMs <- MM$MM$MMs
	MMns <- MM$MM$MMns
	geno.part <- MM$geno.part
	
	ldim <- MM$dim
	random. <- unlist(lapply(ldim, get.attribute, "random"))
	sparse. <- unlist(lapply(ldim, get.attribute, "sparse"))
	spatial. <- unlist(lapply(ldim, get.attribute, "spatial"))
	dim <- unlist(ldim)

	# Nominal dimension
	dim.nom <- obtain.nominal.dimension(cbind(MM$MM$MMs, MM$MM$MMns), dim, random., spatial., weights.na)
	g <- MM$g
	
	random.pos <- create.position.indicator(dim, random.)	
	df.fixed <- sum(dim[!random.])

	# Fixed and random effects
	# bold = rep(0, sum(dim, na.rm = TRUE))
	# Fit the model
	# eta <- cbind.spam(MMs, MMns)%*%bold + offset.na
	# mu <- family$linkinv(eta)
	mustart <- etastart <- NULL
	eval(family$initialize)
	mu <- mustart
	eta <- family$linkfun(mustart)

	# Initialize variance components
	la <- c(1, unlist(MM$init.var))
	# Initialize deviance and psi
	devold <- 1e10
	psi <- la[1]
	if(control$monitoring > 1) {
		cat("Effective dimensions\n")
		cat("-------------------------\n")
		cat(sprintf("%1$3s %2$12s","It.","Deviance"), sep = "")
		cat(sprintf("%12s", names(g)), sep = "")
		cat("\n")
	}
	for (iter in 1:control$maxit) {
		deriv <- family$mu.eta(eta)
		z <- (eta - offset.na) + (y - mu)/deriv
		w <- as.vector(deriv^2/family$variance(mu))
		w <- w*weights.na
		z[!weights.na] <- 0
		mat <- construct.matrices(MMs, MMns, z, w)
		
		if(!genotype.as.random) {
			M <- t((1/mat$XtX.)%*%mat$XtZ.)
			C21_C11_inv_C12 <- M%*%mat$XtZ.
			M_Xty. <- M%*%mat$Xty.
		}
		
		if(control$monitoring) start1 <- proc.time()[3]
		for (it in 1:control$maxit) {
			# Build penalty matrix: block diagonal matrix
			Ginv <- vector(length = sum(dim[random.], na.rm = TRUE))
			for (i in 1:length(g)) {
				Ginv <- Ginv + (1/la[i+1])*g[[i]]
			}
			G <- 1/Ginv
			
			obj <- construct.henderson.matrix(mat, la, Ginv, dim, sparse., random., fixed.matrices = if(genotype.as.random) NULL else list(M = M, C21_C11_inv_C12 = C21_C11_inv_C12, M_Xty. = M_Xty.), as.random = genotype.as.random)
			
			# Fixed and random effects estimation
			chol_K <- try(chol(obj$K), silent = TRUE)
			if(class(chol_K) == "try-error") {
				stop("The design matrix associated to the fixed part of the model is not of full rank. Please check that there are no redundant components in the fixed part of the model.")	
			}

			K_inv <- chol2inv(chol_K)
			
			b.fr <- (1/la[1])*K_inv%*%(mat$Zty. - obj$M_Xty.)
			b.fr <- as.vector(b.fr)
			names(b.fr) <- colnames(MMns)
			
			b.geno <- as.vector(obj$A_inv%*%((1/la[1])*mat$Xty. - obj$C12%*%b.fr))
			names(b.geno) <- attr(MMs, "colnames")

			b <- c(b.geno, b.fr)
			b.random <- b[random.pos]
			
			# Variance components and effective dimensions
			hat_inverses <- compute.hat.diagonal(mat, la, obj$A_inv, K_inv, obj$M, Ginv, dim, sparse., random., as.random = genotype.as.random)
			dZtPZ <- hat_inverses$dZtPZ
			ed <- tau <- vector(mode="list")
			for (i in 1:length(g)) {
				g.inv.d <- (1/la[i+1])*g[[i]]
				ed[[i]] <- sum(dZtPZ*(g.inv.d*G^2))
				ed[[i]] <- ifelse(ed[[i]] == 0, 1e-50, ed[[i]])
				tau[[i]] <- sum(b.random^2*g[[i]])/ed[[i]]
				tau[[i]] <- ifelse(tau[[i]] == 0, 1e-50, tau[[i]])
			}
			#ssr = sum(((z - cbind(MMs, MMns)%*%b)*sqrt(w))^2) #yty. - t(b)%*%(2*u - V%*%b)
			ssr <- sum(((z - MMs %*% b.geno - MMns %*% b.fr)*sqrt(w))^2) #yty. - t(b)%*%(2*u - V%*%b)
			# Compute deviance
			dev <- deviance(obj$C11, chol_K, G, w[w != 0], la[1], ssr, sum(b.random^2*Ginv))
			psinew <- as.numeric((ssr/(nobs - sum(unlist(ed)) - df.fixed)))
			if(family$family == "gaussian" | control$update.psi) {
				psi2 <- psinew
			} else {
				psi2 <- 1
			}
			# New variance components and convergence check
			lanew <- c(psi2, unlist(tau))
			dla <- abs(devold - dev)
			if(control$monitoring > 1) {
				cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
				cat(sprintf("%12.3f", unlist(ed)), sep = "")
				cat('\n')
			}
			if (dla < control$tolerance) break
			la <- lanew
			psi <- psinew
			devold <- dev
		}
		if (control$monitoring) {
			end1 <- proc.time()[3]
			cat("Timings:\nSpATS", (end1-start1), "seconds\n")
		}
		eta.old <- eta
		eta <- cbind.spam(MMs, MMns)%*%b + offset.na
		mu <- family$linkinv(eta)
		# Convergence criterion: linear predictor
		tol <- sum((eta - eta.old)^2)/sum(eta^2)
		if (tol < control$tolerance | (family$family == "gaussian" & family$link== "identity")) break
	}	
	# Inverses
	C11_inv <- as.matrix(obj$A_inv - hat_inverses$C12_inv%*%obj$M)
	C12_inv <- as.matrix(hat_inverses$C12_inv)
	C21_inv <- as.matrix(hat_inverses$C21_inv)
	C22_inv <- as.matrix(hat_inverses$C22_inv)
	
	colnames(C11_inv) <- rownames(C11_inv) <- rownames(C12_inv) <- colnames(C21_inv) <- names(b.geno)
	colnames(C22_inv) <- rownames(C22_inv) <- colnames(C12_inv) <- rownames(C21_inv) <- names(b.fr)
	
	var.comp <- la[-1]
	eff.dim <- unlist(ed)
	names(var.comp) <- names(eff.dim) <- names(g)
	
	# Effective dimension (fixed + random)
	eff.dim <- c(dim[!random.], eff.dim)
	
	attr(dim, "random") <- random.
	attr(dim, "sparse") <- sparse.
	attr(dim, "spatial") <- spatial.
	
	fitted <- rep(NA, nrow(data))
	fitted[!na.ind] <- mu
	
	# Deviance residuals
	dev.residuals <- family$dev.resids(data[,response], fitted, weights)
    s <- attr(dev.residuals, "sign")
    if (is.null(s)) 
        s <- sign(data[,response] - fitted)
    dev.residuals <- sqrt(pmax(dev.residuals, 0))*s
    
	if (control$monitoring) {
		end <- proc.time()[3]
		cat("All process", (end - start), "seconds\n")
	}
	res <- list()
	res$call <- match.call()
	res$data <- cbind(data, weights = weights)
	res$model <- list(response = response, spatial = spatial, geno = list(genotype = genotype, geno.decomp = geno.decomp, as.random = genotype.as.random), fixed = fixed, random = random)
	res$fitted <- fitted
	res$residuals <- dev.residuals
	res$psi <- c(la[1], psi)
	res$var.comp <- var.comp
	res$eff.dim <- eff.dim
	res$dim <- dim
	res$dim.nom <- dim.nom 
	res$nobs <- nobs
	res$deviance <- dev
	res$coeff <- b
	res$niterations <- it
	random.coeff <- rep(FALSE, length(b))
	random.coeff[create.position.indicator(dim, random.)] <- TRUE
	attr(res$coeff, "random") <- random.coeff	
	# Terms
	res$terms <- MM$terms
	# vcov
	res$vcov <- list(C11_inv = C11_inv, C12_inv = C12_inv, C21_inv = C21_inv, C22_inv = C22_inv)
	class(res) <- "SpATS"
	res
}
