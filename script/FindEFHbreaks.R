FindEFHbreaks = function (abund.raster, method = "cumulative", threshold = 0.0513, 
          quantiles = c(0.05, 0.25, 0.5, 0.75), data = NULL) 
{
  quants <- sort(unique(c(0, quantiles, 1)))
  if (method == "percentile") {
    sample <- stats::na.omit(raster::getValues(abund.raster))
    sample[sample <= threshold] <- NA
    breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, 
                              names = FALSE)
    breaks[1] <- 0
    breaks[length(breaks)] <- Inf
  }
  if (method == "cumulative") {
    if (is.null(data)) {
      vals <- stats::na.omit(sort(raster::getValues(abund.raster)))
      vals2 <- cumsum(vals)/sum(vals)
    }
    else {
      vals <- stats::na.omit(sort(raster::extract(abund.raster, 
                                                  data.frame(data$lon, data$lat))))
      vals2 <- cumsum(vals)/sum(vals)
    }
    breaks <- c(0, rep(NA, length(quants) - 2), Inf)
    while (length(unique(stats::na.omit(breaks))) != length(quants)) {
      for (j in 2:(length(quants) - 1)) {
        breaks[j] <- vals[which(vals2 > quants[j])[1]]
      }
      vals <- vals[-length(vals)]
      vals2 <- cumsum(vals)/sum(vals)
    }
  }
  return(breaks)
}