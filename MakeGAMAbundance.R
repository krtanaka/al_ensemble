MakeGAMAbundance = function (model, r.stack, scale.factor = 1, filename = "", land = NULL, 
          mask = NULL) 
{
  if (is.na(filename) | is.null(filename)) {
    filename <- ""
  }
  model.terms <- AutodetectGAMTerms(model)
  if (model$family$family == "ziplss") {
    dterms <- model.terms[[1]]
    pterms <- model.terms[[2]]
    model.terms <- unique(rbind(dterms, pterms))
  }
  if ("offset" %in% model.terms$type) {
    off.name <- model.terms$term[which(model.terms$type == 
                                         "offset")]
    off.val <- ifelse(is.list(model$offset), mean(model$offset[[1]]), 
                      mean(model$offset))
    off.raster <- raster::raster(ext = r.stack@extent, crs = r.stack@crs, 
                                 nrow = r.stack@nrows, ncol = r.stack@ncols, vals = off.val)
    names(off.raster) <- off.name
    r.stack <- raster::stack(list(r.stack, off.raster))
  }
  gam.factors <- model.terms$term[model.terms$type == "factor"]
  gam.factors2 <- list()
  if (length(gam.factors) > 0) {
    for (t in 1:length(gam.factors)) {
      range <- raster::subset(x = r.stack, subset = which(names(r.stack) == 
                                                            gam.factors[t]))
      gam.factors2[[t]] <- sort(unique(stats::na.omit(raster::getValues(range))))
    }
    names(gam.factors2) <- gam.factors
  }
  else {
    gam.factors2 <- NULL
  }
  pred.type <- ifelse(model$family$link == "cloglog", "link", 
                      "response")[1]
  if (model$family$family == "ziplss" & is.null(gam.factors2)) {
    r.vals <- as.data.frame(raster::getValues(r.stack))
    if ("offset" %in% model.terms$type) {
      r.vals <- cbind(r.vals, data.frame(off.val))
      names(r.vals[ncol(r.vals)]) <- off.name
    }
    pred.vals <- mgcv::predict.gam(model, newdata = r.vals, 
                                   type = "response")
    predict.raster <- raster::setValues(x = raster::raster(r.stack), 
                                        values = pred.vals)
  }
  else {
    predict.raster <- raster::predict(r.stack, model, factors = gam.factors2, 
                                      progress = "text", overwrite = TRUE, type = pred.type, 
                                      newdata.guaranteed = TRUE)
  }
  if (model$family$link[1] == "cloglog") {
    predict.raster <- exp(predict.raster)
  }
  predict.raster <- predict.raster * scale.factor
  if (filename != "") {
    raster::writeRaster(x = predict.raster, filename = filename, 
                        overwrite = TRUE)
  }
  if (is.null(land) == F) {
    predict.raster <- raster::mask(predict.raster, land, 
                                   inverse = TRUE, overwrite = TRUE, filename = filename)
  }
  if (is.null(mask) == F) {
    predict.raster <- raster::mask(predict.raster, mask, 
                                   overwrite = TRUE, filename = filename)
  }
  return(predict.raster)
}