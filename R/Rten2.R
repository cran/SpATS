Rten2 <-
function(X1,X2) {
	one.1 <- matrix(1,1,ncol(X1))
	one.2 <- matrix(1,1,ncol(X2))
	kronecker(X1,one.2)*kronecker(one.1,X2)
}
