FitMaxnet <- function(data, species, vars, reduce = FALSE, regmult = 1, facs = NULL) {
  
  presence.vec <- as.integer(data[, species] > 0)
  maxnet.data <- data[, c(vars, facs)]
  
  if (length(facs) > 0) {
    for (f in facs) {
      maxnet.data[, f] <- as.factor(maxnet.data[, f])
    }
  }
  
  drops <- NULL
  for (i in 1:ncol(maxnet.data)) {
    drops <- c(drops, which(is.na(maxnet.data[, i])))
  }
  
  if (length(drops) > 0) {
    unique_drops <- unique(drops)
    presence.vec <- presence.vec[-unique_drops]
    maxnet.data <- maxnet.data[-unique_drops, ]
  }
  
  maxnet.model <- maxnet::maxnet(p = presence.vec, data = maxnet.data, regmult = regmult)
  
  if (reduce) {
    m.coefs <- MaxnetCoefs(maxnet.model)
    badvars <- names(m.coefs)[which(m.coefs == 0)]
    
    if (length(badvars) > 0) {
      badcols <- which(names(maxnet.data) %in% badvars)
      maxnet.data2 <- maxnet.data[, -badcols]
      maxnet.model <- maxnet::maxnet(p = presence.vec, data = maxnet.data2, regmult = regmult)
    }
  }
  
  return(maxnet.model)
}
