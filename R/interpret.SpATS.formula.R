interpret.SpATS.formula <-
function(formula) {
    env <- environment(formula) 
    if(inherits(formula, "character"))          
        formula <- as.formula(formula)
    tf <- terms.formula(formula, specials = c("SAP", "PSANOVA"))
    terms <- attr(tf, "term.labels")
    nt <- length(terms)
    if(nt != 1)
    	stop("Error in the specification of the spatial effect: only a sigle bidimensional function is allowed")
    
    res <- eval(parse(text = terms[1]), envir = env)
    res
}
