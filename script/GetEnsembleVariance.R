GetEnsembleVariance = function (model.weights, variance.list, abund.list, ensemble.abund) 
{
  keepers <- which(model.weights > 0)
  weights2 <- model.weights[keepers]
  variance.list2 <- list()
  abund.list2 <- list()
  for (k in 1:length(keepers)) {
    variance.list2[[k]] <- variance.list[[keepers[k]]]
    abund.list2[[k]] <- abund.list[[keepers[k]]]
  }
  data.spots <- which(is.na(raster::getValues(abund.list2[[1]])) == 
                        F)
  if (length(abund.list2) > 1) {
    for (r in 2:length(abund.list2)) {
      spots <- which(is.na(raster::getValues(abund.list2[[r]])) == 
                       F)
      data.spots <- data.spots[data.spots %in% spots]
    }
  }
  e.abund <- raster::extract(ensemble.abund, data.spots)
  dat <- matrix(nrow = length(data.spots), ncol = length(keepers))
  for (m in 1:length(keepers)) {
    a.dat <- raster::extract(abund.list2[[m]], data.spots)
    v.dat <- raster::extract(variance.list2[[m]], data.spots)
    dat[, m] <- weights2[m] * sqrt(v.dat + (a.dat - e.abund)^2)
  }
  std.error <- apply(dat, MARGIN = 1, FUN = sum)
  out.raster <- raster::raster(abund.list2[[1]])
  val.vec <- raster::getValues(out.raster)
  val.vec[data.spots] <- std.error
  out.raster <- raster::setValues(out.raster, val.vec)
  return(out.raster)
}