deviance <-
function(A, K, G, w, sigma2, ssr, edf) {
	log_det_A <- sum(log(A))
	log_det_K <- 2*sum(log(diag(K)))
	log_det_G <- sum(log(G))
	deviance <- log_det_A + log_det_K + log_det_G + sum(log(sigma2*1/w)) + ssr/sigma2 + edf
	deviance	
}
