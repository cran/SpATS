variogram.SpATS <-
function(x, ...) {
	xlab <- x$terms$spatial$terms.formula$x.coord
	ylab <- x$terms$spatial$terms.formula$y.coord
	x.coord <- x$data[,xlab]
	y.coord <- x$data[,ylab]
	residuals <- x$residuals
	
	columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
	rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
	
	xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
	setkeyv(xy.coord, c("rows", "columns"))
	df <- data.table(columns = x.coord, rows = y.coord, residuals = residuals)
	setkeyv(df, c("rows", "columns"))
	df <- df[xy.coord]
	df <- df[order(df$columns, df$rows),]
	
	resdiff <- c(outer(df$residuals, df$residuals, function(x,y) 0.5*(x-y)^2))
	coldiff <- c(outer(df$columns, df$columns, function(x,y) abs(x-y)))
	coldiff.u <- unique(coldiff)
	rowdiff <- c(outer(df$rows, df$rows, function(x,y) abs(x-y)))
	rowdiff.u <- unique(rowdiff)

	subsets <- split(resdiff, f = list(coldiff, rowdiff))
	value <- sapply(subsets, mean, na.rm = TRUE)
	length <-  sapply(subsets, function(x) sum(!is.na(x)))
	length[-1] <- length[-1]/2

	res <- list(data = data.frame(value = value, length = length), col.displacement = coldiff.u, row.displacement = rowdiff.u)
	class(res) <- "variogram.SpATS"
	res

}
