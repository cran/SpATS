construct.capital.lambda <-
function(g) {
	length.eq <- all(sapply(g, function(x) {
		 	diff(range(unlist(lapply(x, length)))) < .Machine$double.eps ^ 0.5
		}))
	if(length.eq) {
		l <- length(g)
		if(l == 1) {
			if(length(g[[1]]) == 1) {
				res <- g
			 } else {
			 	res <- do.call("c", lapply(g, function(x) x))
			 }
		} else {
			dim <- sapply(g, function(x) {
				if(is.list(x))
		 			unlist(lapply(x, length))[1]
		 		else
		 			length(x)
			})		
			end <- cumsum(dim)
			init <- end - dim + 1
				
			res <- do.call("c", lapply(1:length(g), function(x, g, init, end, dim) {
				temp <- g[[x]]
				if(is.list(temp)) {
					lapply(temp, function(y, x, dim) {
						aux <- rep(0, l = sum(dim))
						aux[init[x]:end[x]] <- y
						aux
					}, x = x, dim = dim)
				} else {
					aux <- rep(0, l = sum(dim))
					aux[init[x]:end[x]] <- temp
					list(aux)
				}
			}, g = g, init = init, end = end, dim = dim))	
		}
	} else {
		stop("Error in construct.capital.lambda")
	}	
	res
}
