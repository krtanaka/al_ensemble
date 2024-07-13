MakeEnsemble = function (ensemble.list = NA, rmse = NA, names.vec = NA, minimum = NA, 
          model.types = NULL) 
{
  if (any(is.na(rmse))) {
    if (is.na(names.vec)) {
      names.vec <- names(ensemble.list)
    }
    rmse <- vector(length = length(ensemble.list))
    for (m in 1:length(ensemble.list)) {
      ensemble.dat <- stats::na.omit(ensemble.list[[m]])
      if (nrow(ensemble.dat) > 0) {
        rmse[m] <- sqrt(sum((ensemble.dat$abund - ensemble.dat$pred)^2)/nrow(ensemble.dat))
      }
      else {
        rmse[m] <- NA
      }
    }
  }
  else {
    if (is.na(names.vec)) {
      names.vec <- names(rmse)
    }
  }
  rmse2 <- rmse^2
  weights <- (1/(rmse2))/sum(1/(rmse2), na.rm = T)
  names(weights) <- names.vec
  weights[is.na(weights)] <- 0
  if (length(model.types) > 0) {
    for (m in unique(model.types)) {
      index <- which(model.types == m)
      weights2 <- weights[index]
      keep <- index[which.max(weights2)]
      drop <- index[index != keep]
      weights[drop] <- 0
    }
    weights <- weights/sum(weights)
  }
  if (is.na(minimum) == F) {
    included <- which(weights >= minimum)
    set.to.zero <- which(weights < minimum)
    weights2 <- weights[included]
    weights2 <- weights2/sum(weights2)
    weights[set.to.zero] <- 0
    weights[included] <- weights2
  }
  return(weights)
}