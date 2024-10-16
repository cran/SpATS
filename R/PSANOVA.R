PSANOVA <-
function(..., nseg = c(10, 10), pord = c(2,2), degree = c(3,3), nest.div = c(1,1), center = FALSE) {   
    args <- match.call()
	vars <- as.list(substitute(list(...)))[-1]
    if(length(vars) != 2) {
        stop("Error in the specification of the spatial effect: both spatial coordinates must be indicated")
    }
  	
  	if(any(pord != 2)) {
  		stop("Error in the specification of the spatial effect: Only second order penalties are allowed in PSANOVA")
    }

    nseg <- if(length(nseg) == 1) {
        rep(nseg, 2) 
    } else {
        nseg
    }

    pord <- if(length(pord) == 1) {
        rep(pord, 2) 
    } else {
        pord
    }

    nest.div <- if(length(nest.div) == 1) {
        rep(nest.div, 2)
    } else { 
        nest.div
    }

    if(any(nseg%%nest.div != 0)) {
        stop("The number of segments, <nseg>, must be a multiple of the divisor of the number of segments, <nest.div>.")
    }

  		
  	x.coord = vars[[1]]
  	y.coord = vars[[2]]
  	
  	res <- list()
  	res$x.coord <- deparse(x.coord, backtick = TRUE, width.cutoff = 500)
  	res$y.coord <- deparse(y.coord, backtick = TRUE, width.cutoff = 500)
  	res$nseg <- nseg
  	res$pord <- pord
  	res$degree <- if(length(degree) == 1) rep(degree, 2) else degree
  	res$nest.div <- nest.div
    res$center <- center
  	res$type <- "PSANOVA"
  	
  	res
}
