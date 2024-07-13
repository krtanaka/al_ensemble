MakeVarianceRasters = function (model.list, raster.stack, model.type, scale.factor = 1, efh.break = NA) {
  
  data.spots <- which(is.na(raster::getValues(raster.stack[[1]])) == FALSE)
  
  for (r in 2:raster::nlayers(raster.stack)) {
    spots <- which(is.na(raster::getValues(raster.stack[[r]])) == FALSE)
    data.spots <- data.spots[data.spots %in% spots]
  }
  
  data <- extract(raster.stack, data.spots)
  
  list.index <- 1
  model.list2 <- list()
  
  for (m in 1:length(model.list)) {
    if (is.list(model.list[[m]])) {
      model.list2[[list.index]] <- model.list[[m]]
      list.index <- list.index + 1
    }
  }
  
  out.data <- matrix(nrow = nrow(data), ncol = length(model.list2))
  
  if (model.type != "maxnet") {
    terms <- AutodetectGAMTerms(model.list[[1]], hgam = "d")
    has.offset <- "offset" %in% terms$type
    facs <- terms$term[terms$type == "factor"]
    offset.name <- terms$term[terms$type == "offset"]
    
    if (model.type == "hgam" & has.offset) {
      data <- data.frame(data, mean(model.list2[[1]]$offset[[1]]))
      pterms <- AutodetectGAMTerms(model.list[[1]], hgam = "p")
      facs <- unique(c(facs, pterms$term[pterms$type == "factor"]))
    }
    
    if (model.type %in% c("cloglog", "gam")) {
      data <- data.frame(data, mean(model.list2[[1]]$offset))
    }
    
    names(data)[ncol(data)] <- offset.name
    
    for (f in 1:length(facs)) {
      data[, facs[f]] <- as.factor(data[, facs[f]])
    }
  }
  
  pb <- utils::txtProgressBar(min = 0, max = length(model.list2), style = 3)
  
  for (m in 1:length(model.list2)) {
    if (model.type == "maxnet") {
      out.data[, m] <- exp(stats::predict(model.list2[[m]], newdata = data, type = "link") + model.list2[[m]]$entropy)
    }
    if (model.type == "cloglog") {
      out.data[, m] <- exp(mgcv::predict.gam(object = model.list2[[m]], newdata = data, type = "link"))
    }
    if (model.type == "hgam") {
      out.data[, m] <- mgcv::predict.gam(object = model.list2[[m]], newdata = data, type = "response")
    }
    if (model.type == "gam") {
      out.data[, m] <- mgcv::predict.gam(object = model.list2[[m]], newdata = data, type = "response")
    }
    utils::setTxtProgressBar(pb, m)
  }
  
  close(pb)
  
  raster.template <- raster::raster(raster.stack)
  variances <- apply(X = out.data * scale.factor, MARGIN = 1, FUN = stats::var)
  var.vec <- rep(NA, times = raster::ncell(raster.stack))
  var.vec[data.spots] <- variances
  var.raster <- raster::setValues(raster.template, values = var.vec)
  
  if (!is.na(efh.break)) {
    percents <- apply(X = out.data > efh.break, MARGIN = 1, FUN = sum) / ncol(out.data)
    per.vec <- rep(NA, times = raster::ncell(raster.stack))
    per.vec[data.spots] <- percents
    per.raster <- raster::setValues(raster.template, values = per.vec)
    return(list(var.raster, per.raster))
  } else {
    return(var.raster)
  }
  
}