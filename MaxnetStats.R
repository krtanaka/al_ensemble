MaxnetStats = function (model, maxnet2d = NULL, regmult = 1, data, species) 
{
  covar.names <- names(model$samplemeans)
  if (length(maxnet2d) > 0) {
    coef.vec2 <- vector(length = length(maxnet2d))
    vars1d <- covar.names[covar.names %in% unlist(maxnet2d) == 
                            F]
    maxnet2d.name <- unlist(lapply(X = maxnet2d, FUN = function(x) {
      paste(x, collapse = "*")
    }))
  }
  else {
    vars1d <- covar.names
  }
  min.names <- names(model$varmin)
  facs <- vars1d[vars1d %in% min.names == F]
  vars1d <- vars1d[vars1d %in% min.names]
  pb <- utils::txtProgressBar(min = 0, max = length(c(maxnet2d, 
                                                      vars1d, facs)), style = 3)
  pb.i <- 1
  if (length(maxnet2d) > 0) {
    dev.vec2d <- rep(NA, length = length(maxnet2d))
    for (i in 1:length(maxnet2d)) {
      d.covars <- covar.names[covar.names %in% maxnet2d[[i]] == 
                                F]
      try(test.model <- FitMaxnet(data = data, species = species, 
                                  vars = d.covars, facs = facs, regmult = regmult))
      if (exists("test.model")) {
        dev.vec2d[i] <- test.model$dev.ratio[length(test.model$dev.ratio)]
        rm(test.model)
        utils::setTxtProgressBar(pb, pb.i)
        pb.i <- pb.i + 1
      }
      else {
        close(pb)
        break
      }
    }
    names(dev.vec2d) <- maxnet2d.name
  }
  else {
    dev.vec2d <- NULL
    maxnet2d.name <- NULL
  }
  dev.vec <- rep(NA, length = length(vars1d))
  for (i in 1:length(vars1d)) {
    try(test.model <- FitMaxnet(data = data, species = species, 
                                vars = c(unlist(maxnet2d), vars1d[-i]), facs = facs, 
                                regmult = regmult))
    if (exists("test.model")) {
      dev.vec[i] <- test.model$dev.ratio[length(test.model$dev.ratio)]
      rm(test.model)
      utils::setTxtProgressBar(pb, pb.i)
      pb.i <- pb.i + 1
    }
    else {
      close(pb)
      break
    }
  }
  names(dev.vec) <- vars1d
  fac.vec <- rep(NA, length = length(facs))
  for (i in 1:length(facs)) {
    try(test.model <- FitMaxnet(data = data, species = species, 
                                vars = c(unlist(maxnet2d), vars1d), facs = facs[-i], 
                                regmult = regmult))
    if (exists("test.model")) {
      fac.vec[i] <- test.model$dev.ratio[length(test.model$dev.ratio)]
      rm(test.model)
      utils::setTxtProgressBar(pb, pb.i)
      pb.i <- pb.i + 1
    }
    else {
      close(pb)
      break
    }
  }
  close(pb)
  names(fac.vec) <- facs
  out.dev.vec <- c(dev.vec2d, dev.vec, fac.vec)
  if (sum(is.na(out.dev.vec)) == 0) {
    dev.lost <- 1 - out.dev.vec/model$dev.ratio[length(model$dev.ratio)]
    dev.exp <- dev.lost/sum(dev.lost, na.rm = T) * 100
    if (min(dev.exp) < 0) {
      dev.exp <- dev.exp - min(dev.exp)
      dev.exp <- dev.exp/sum(dev.exp) * 100
    }
  }
  else {
    dev.exp <- rep(NA, times = length(c(maxnet2d, vars1d, 
                                        facs)))
  }
  names(dev.exp) <- c(maxnet2d.name, vars1d, facs)
  return(dev.exp)
}