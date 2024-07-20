MakeEnsembleAbundance = function (model.weights, abund.list, filename = "") 
{
  good.abund <- which(is.na(abund.list) == F)
  new.abund <- raster::raster(abund.list[[good.abund[1]]])
  new.abund <- raster::setValues(new.abund, values = 0)
  for (i in 1:length(model.weights)) {
    if (model.weights[i] > 0) {
      new.abund <- new.abund + abund.list[[i]] * model.weights[i]
    }
  }
  if (filename != "") {
    raster::writeRaster(x = new.abund, filename = filename, 
                        overwrite = T)
  }
  return(new.abund)
}