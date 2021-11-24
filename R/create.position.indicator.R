create.position.indicator <-
function(dim, what) {
    end.tot <- cumsum(dim)
    init.tot <- end.tot - dim + 1
    
    init <- init.tot[what]
    end <- end.tot[what]
        
    if((length(init) == 0 | is.null(init)) & (length(end) == 0 | is.null(end))) {
        res <- NULL
    } else if(length(init) != 0 & length(end) != 0) {
        res <- unlist(sapply(1:length(init), function(i, init, end) {
            res <- init[i]:end[i]
            res
        }, init = init, end = end))
    } else {
        stop("Error in create.position.indicator")
    }    
    as.vector(res)
}
