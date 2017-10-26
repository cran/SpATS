plot.variogram.SpATS <-
function(x, min.length = 30, ...) {
	values <- matrix(replace(x$data$value, x$data$length < min.length, NA), ncol = length(x$col.displacement), nrow = length(x$row.displacement), byrow = TRUE)
	plot3Drgl::persp3Drgl(x$row.displacement, x$col.displacement, values, xlab = "Row displacement", ylab = "Col displacement", zlab = "", ticktype = "detailed", col = topo.colors(100))
}
