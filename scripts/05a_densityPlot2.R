densityPlot2<-function (dat, sampGroups = NULL, main = "", xlab = "Beta", pal = brewer.pal(8, 
                                                                             "Dark2"), xlim, ylim, add = TRUE, legend = TRUE, ...) 
{
  if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
    b <- getBeta(dat)
  }
  else if (is(dat, "matrix")) {
    b <- dat
  }
  else {
    stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet' or ", 
         "matrix.")
  }
  d <- apply(b, 2, function(x) density(as.vector(x), na.rm = TRUE))
  if (missing(ylim)) 
    ylim <- range(sapply(d, function(i) range(i$y)))
  if (missing(xlim)) 
    xlim <- range(sapply(d, function(i) range(i$x)))
  if (is.null(sampGroups)) {
    sampGroups <- rep(1, ncol(b))
  }
  else if (length(sampGroups) == 1) {
    sampGroups <- rep(sampGroups, ncol(b))
  }
  sampGroups <- as.factor(sampGroups)
  if (add) {
    plot(x = 0, type = "n", ylim = ylim, xlim = xlim, ylab = "Density", 
         xlab = xlab, main = main, ...)
    abline(h = 0, col = "grey80")
  }
  for (i in seq_along(d)) {
    lines(d[[i]], col = pal[sampGroups[i]])
  }
  if (legend & length(levels(sampGroups)) > 1) {
    legend("top", legend = levels(sampGroups), text.col = pal)
  }
}