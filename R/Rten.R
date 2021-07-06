Rten <-
function(X) {
	one <- matrix(1, 1, ncol(X))
	kronecker(X,one)*kronecker(one,X)
}
