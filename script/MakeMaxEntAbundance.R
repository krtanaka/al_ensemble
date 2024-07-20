MakeMaxEntAbundance <- function(model, maxent.stack, scale.fac = 1, land = NULL, mask = NULL, 
                                type = "cloglog", clamp = FALSE, filename = "") {
  if (is.null(filename) || is.na(filename)) {
    filename <- ""
  }
  
  type <- tolower(type)
  
  dat <- raster::getValues(maxent.stack)
  na.spots <- which(apply(X = dat, MARGIN = 1, FUN = function(x) {
    return(any(is.na(x)))
  }))
  dat.spots <- which(seq(1:nrow(dat)) %in% na.spots == FALSE)
  
  if (type == "cloglog") {
    preds <- stats::predict(model, dat[dat.spots, ], type = "link")
    preds2 <- exp(preds + model$ent) * scale.fac
  } else if (type == "maxnet") {
    preds <- stats::predict(model, dat[dat.spots, ], type = "cloglog")
    preds2 <- preds
  } else {
    stop("Unknown type: ", type)
  }
  
  new.vals <- vector(length = nrow(dat))
  new.vals[na.spots] <- NA
  new.vals[dat.spots] <- preds2
  
  habitat.prediction <- raster::setValues(x = raster::raster(maxent.stack), values = new.vals)
  
  if (!is.null(land)) {
    habitat.prediction <- raster::extend(x = habitat.prediction, y = land)
  }
  
  habitat.prediction@crs <- maxent.stack@crs
  
  if (filename != "") {
    raster::writeRaster(x = habitat.prediction, filename = filename, overwrite = TRUE)
  }
  
  if (!is.null(land)) {
    habitat.prediction <- raster::mask(habitat.prediction, land, inverse = TRUE, overwrite = TRUE, filename = filename)
  }
  
  if (!is.null(mask)) {
    habitat.prediction <- raster::mask(habitat.prediction, mask, overwrite = TRUE, filename = filename)
  }
  
  return(habitat.prediction)
}
