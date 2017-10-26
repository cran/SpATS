plot.SpATS <-
function(x, all.in.one = TRUE, main = NULL, annotated = FALSE, depict.missing = FALSE, ...) {
    xlab <- x$terms$spatial$terms.formula$x.coord
    ylab <- x$terms$spatial$terms.formula$y.coord
    x.coord <- x$data[,xlab]
    y.coord <- x$data[,ylab]
    response <- x$data[,x$model$response]
    residuals <- x$residuals
    fitted <- x$fitted
    geno.model.matrix <- construct.genotype.prediction.matrix(x, x$data)
    geno.coeff <- x$coeff[1:ncol(geno.model.matrix)]
    geno.pred <- as.vector(geno.model.matrix%*%geno.coeff)
    
    # NAs
    residuals[x$data$weights == 0] <- NA
    fitted[x$data$weights == 0] <- NA
    geno.pred[x$data$weights == 0] <- NA
    
    if(is.null(main)) main = paste("Trait: ", x$model$response, sep = "")
    
    columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
    rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
    
    setNumericRounding(2)
    
    xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
    setkeyv(xy.coord, c("rows", "columns"))
    ONE <- rep(1, length(x.coord))    
    df <- data.table(columns = x.coord, rows = y.coord, response = response, fitted = fitted, residuals = residuals, geno.pred = geno.pred, weights = x$data$weights, ONE = ONE)
    setkeyv(df, c("rows", "columns"))
    df <- df[xy.coord]
    df <- df[order(df$columns, df$rows),]
    if(depict.missing)
        df$ONE[is.na(df$ONE)] <- 1
    else
        df$ONE[df$weights == 0 | is.na(df$weights)] <- NA
        
        
    p1 <- if(length(columns) > 100) 1 else 100%/%length(columns) + 1
    p2 <- if(length(rows) > 100) 1 else 100%/%length(rows) + 1

    fit.spatial.trend <- obtain.spatialtrend(x, grid = c(length(columns)*p1, length(rows)*p2))
    Mf = kronecker(matrix(df$ONE, ncol = length(columns), nrow = length(rows)), matrix(1, p2, p1))
    
    colors = topo.colors(100)
    
    main.legends <- c('Raw data', 'Fitted data', 'Residuals', 'Fitted Spatial Trend', ifelse(x$model$geno$as.random, "Genotypic BLUPs", "Genotypic BLUEs"), 'Histogram')
    if(all.in.one) {
        op <- par(mfrow = c(2,3), oma = c(ifelse(annotated, 12, 2), 1, 3, 2), mar = c(2.5, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))                
    } else {
        if(!is.null(main))
            main.legends <- rep(main, length(main.legends))
    }

    range <- range(c(response, fitted), na.rm = TRUE)
    fields::image.plot(columns, rows, t(matrix(df$response, ncol = length(columns), nrow = length(rows))), main = main.legends[1], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
    if(!all.in.one)
        readline("Press return for next page....")
    fields::image.plot(columns, rows, t(matrix(df$fitted, ncol = length(columns), nrow = length(rows))), main = main.legends[2], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
    if(!all.in.one)
        readline("Press return for next page....")
    fields::image.plot(columns, rows, t(matrix(df$residuals, ncol = length(columns), nrow = length(rows))), main = main.legends[3], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
    if(!all.in.one)
        readline("Press return for next page....")
    fields::image.plot(fit.spatial.trend$col.p, fit.spatial.trend$row.p, t(fit.spatial.trend$fit*Mf), main = main.legends[4], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
    if(!all.in.one)
        readline("Press return for next page....")
    fields::image.plot(columns, rows, t(matrix(df$geno.pred, ncol = length(columns), nrow = length(rows))), main = main.legends[5], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
    if(!all.in.one)
        readline("Press return for next page....")
    suppressWarnings(hist(geno.coeff, main = main.legends[6], xlab = main.legends[5], ...))        
        
    if(all.in.one) {
        title("")
        mtext(main, cex = 1.5, outer = TRUE, side = 3)
        if(annotated) {
            genotype.lab <- ifelse(x$model$geno$as.random, "Genotypes (as random):", "Genotypes (as fixed):")
            text <- paste('\nSpatial analysis of trials with splines \n\n', 
                  sprintf('%-25s %-10s', 'Response:', x$model$response), '\n',
                  sprintf('%-25s %-10s', genotype.lab, x$model$geno$genotype), '\n',
                  sprintf('%-25s %-10s', 'Spatial:', Reduce(paste, deparse(x$model$spatial, width.cutoff = 200L))),
                  if(!is.null(x$model$fixed)) paste('\n', sprintf('%-25s %-10s', 'Fixed:', Reduce(paste, deparse(x$model$fixed, width.cutoff = 200L))), sep = ' ') else NULL,
                  if(!is.null(x$model$random)) paste('\n', sprintf('%-25s %-10s', 'Random:', Reduce(paste, deparse(x$model$random, width.cutoff = 200L))), sep = ' ') else NULL, sep = ' ')
            mtext(text, cex = 0.9, outer = TRUE, side = 1, line = 10, adj = 0, family = "mono")
        }
        par(op)
    }
    invisible(df)
}