plot.variogram.SpATS <-
function(x, min.length = 30, ...) {
	values <- matrix(replace(x$data$value, x$data$length < min.length, NA), ncol = length(x$col.displacement), nrow = length(x$row.displacement), byrow = TRUE)
	#plot3Drgl::persp3Drgl(x$row.displacement, x$col.displacement, values, xlab = "Row displacement", ylab = "Col displacement", zlab = "", ticktype = "detailed", col = topo.colors(100))

	# plot3Drgl no longer works!! 
	nrz <- nrow(values)
	ncz <- ncol(values)
	# Create a function interpolating colors in the range of specified colors
	color <- topo.colors(100)
	# Compute the z-value at the facet centres
	zfacet <- values[-1, -1] + values[-1, -ncz] + values[-nrz, -1] + values[-nrz, -ncz]
	facetcol <- cut(zfacet, 100)

	persp(x$row.displacement, x$col.displacement, values, xlab = "Row displacement", ylab = "Col displacement", zlab = "", ticktype = "detailed", theta = 30, phi = 30, expand = 0.5, col = color[facetcol])
}
