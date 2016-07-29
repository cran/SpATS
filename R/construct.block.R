construct.block <-
function(A1,A2,A3,A4) {
	block <- rbind(cbind(A1,A2), cbind(A3,A4))
	return(block)
}
