controlSpATS <-
function(
	maxit = 200,
	tolerance = 1e-3,
	monitoring = 2,
	update.psi = FALSE,
	update.psi.gauss = TRUE)
	list(maxit = maxit, tolerance = tolerance, monitoring = monitoring, update.psi = update.psi, update.psi.gauss = update.psi.gauss)
