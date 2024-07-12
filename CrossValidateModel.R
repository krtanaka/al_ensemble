CrossValidateModel <- function(
    model, 
    regmult = 1, 
    model.type, 
    scale.preds = FALSE, 
    data, 
    key = NA, 
    species = NA, 
    folds = 10, 
    group = "random"
) {
  
  if (model.type != "maxnet") {
    species <- ifelse(
      model$family$family == "ziplss", 
      as.character(stats::formula(model)[[1]])[[2]], 
      as.character(stats::formula(model))[[2]]
    )
  }
  
  if (!model.type %in% c("maxnet", "cloglog", "hgam", "gam")) {
    stop("Model type not recognized")
  }
  
  model.list <- list()
  scale.factor <- 1
  
  if (group == "random") {
    n.tot <- nrow(data)
    n.pres <- sum(data[, species] > 0)
    pres <- which(data[, species] > 0)
    n.abs <- sum(data[, species] == 0)
    abs <- which(data[, species] == 0)
    pres.base <- floor(n.pres / folds)
    pres.add <- n.pres %% folds
    abs.base <- floor(n.abs / folds)
    abs.add <- n.abs %% folds
    pres.scheme <- c(rep(pres.base + 1, times = pres.add), rep(pres.base, times = folds - pres.add))
    abs.scheme <- c(rep(abs.base, times = folds - abs.add), rep(abs.base + 1, times = abs.add))
    pres.randos <- sample(1:n.pres, size = n.pres, replace = FALSE)
    abs.randos <- sample(1:n.abs, size = n.abs, replace = FALSE)
    pres.group <- rep(LETTERS[1], times = pres.scheme[1])
    abs.group <- rep(LETTERS[1], times = abs.scheme[1])
    
    for (i in 2:folds) {
      pres.group <- c(pres.group, rep(LETTERS[i], times = pres.scheme[i]))
      abs.group <- c(abs.group, rep(LETTERS[i], times = abs.scheme[i]))
    }
    
    pres.data <- data.frame(data[pres[pres.randos], ], group = pres.group)
    abs.data <- data.frame(data[abs[abs.randos], ], group = abs.group)
    data <- rbind(pres.data, abs.data)
    group <- "group"
  }
  
  n.folds <- length(unique(data[, group]))
  fold.vec <- sort(unique(data[, group]))
  fold.table <- table(data[, group])
  start.vec <- cumsum(c(1, fold.table[1:(length(fold.table) - 1)]))
  names(start.vec) <- names(fold.table)
  end.vec <- cumsum(fold.table)
  out.names <- NULL
  
  if (!is.na(key)) {
    out.names <- key
  }
  
  if (sum(c("lon", "lat") %in% names(data)) == 2) {
    out.names <- c(out.names, "lon", "lat")
  }
  
  out.names <- c(out.names, group, "abund", "pred", "prob", "cvpred", "cvprob", "error", "cverror")
  error.data <- as.data.frame(matrix(data = NA, nrow = nrow(data), ncol = length(out.names)))
  colnames(error.data) <- out.names
  
  if (scale.preds) {
    scale.factors <- rep(NA, length = n.folds)
  }
  
  pb <- utils::txtProgressBar(min = 0, max = n.folds, style = 3)
  
  for (i in 1:n.folds) {
    train.data <- data[data[, group] != fold.vec[i], ]
    test.data <- data[data[, group] == fold.vec[i], ]
    error.data[start.vec[i]:end.vec[i], group] <- fold.vec[i]
    
    if (!is.na(key)) {
      error.data[start.vec[i]:end.vec[i], 1] <- test.data[, key]
    }
    
    if ("lon" %in% out.names) {
      error.data$lon[start.vec[i]:end.vec[i]] <- test.data$lon
      error.data$lat[start.vec[i]:end.vec[i]] <- test.data$lat
    }
    
    error.data[start.vec[i]:end.vec[i], group] <- fold.vec[i]
    error.data$abund[start.vec[i]:end.vec[i]] <- test.data[, species]
    
    if (model.type == "maxnet") {
      preds <- exp(predict(object = model, newdata = test.data, response = "link") + model$entropy)
      probs <- predict(object = model, newdata = test.data, type = "cloglog")
      vars0 <- names(model$samplemeans)
      facs <- vars0[!vars0 %in% names(model$varmax)]
      
      try(cv.model <- FitMaxnet(data = train.data, species = species, vars = names(model$varmax), facs = facs, regmult = regmult))
      
      if (exists("cv.model")) {
        cvpreds <- exp(predict(object = cv.model, newdata = test.data, response = "link") + cv.model$entropy)
        cvprobs <- predict(object = cv.model, newdata = test.data, type = "cloglog")
      } else {
        cvpreds <- rep(NA, times = nrow(test.data))
        cvprobs <- rep(NA, times = nrow(test.data))
        cv.model <- NA
      }
    }
    
    if (model.type == "cloglog") {
      preds <- exp(mgcv::predict.gam(object = model, newdata = test.data, type = "link"))
      probs <- mgcv::predict.gam(object = model, newdata = test.data, type = "response")
      
      try(cv.model <- FitGAM(data = train.data, reduce = FALSE, family.gam = "binomial", select = FALSE, link.fx = "cloglog", gam.formula = stats::formula(model), verbose = FALSE))
      
      if (exists("cv.model")) {
        cvpreds <- exp(mgcv::predict.gam(object = cv.model, newdata = test.data, type = "link"))
        cvprobs <- mgcv::predict.gam(object = cv.model, newdata = test.data, type = "response")
      } else {
        cvpreds <- rep(NA, times = nrow(test.data))
        cvprobs <- rep(NA, times = nrow(test.data))
        cv.model <- NA
      }
    }
    
    if (model.type == "hgam") {
      preds <- mgcv::predict.gam(object = model, newdata = test.data, type = "response")
      probs <- 1 - exp(-exp(mgcv::predict.gam(model, newdata = test.data)[, 2]))
      
      try(cv.model <- FitHurdleGAM(density.formula = stats::formula(model)[[1]], prob.formula = stats::formula(model)[[2]], data = train.data, reduce = FALSE, verbose = FALSE, select = FALSE))
      
      if (exists("cv.model")) {
        cvpreds <- mgcv::predict.gam(object = cv.model, newdata = test.data, type = "response")
        cvprobs <- 1 - exp(-exp(mgcv::predict.gam(cv.model, newdata = test.data)[, 2]))
      } else {
        cvpreds <- rep(NA, times = nrow(test.data))
        cvprobs <- rep(NA, times = nrow(test.data))
        cv.model <- NA
      }
    }
    
    if (model.type == "gam") {
      if (strsplit(model$family$family, split = "[()]")[[1]][1] == "Negative Binomial") {
        gamfam <- "nb"
        theta <- as.numeric(strsplit(model$family[[1]], split = "[()]")[[1]][2])
        probs <- 1 - stats::dnbinom(0, mu = mgcv::predict.gam(model, newdata = test.data, type = "response"), size = theta)
      } else {
        gamfam <- model$family$family
        probs <- 1 - stats::dpois(0, mgcv::predict.gam(object = model, newdata = test.data, type = "response"))
      }
      
      preds <- mgcv::predict.gam(object = model, newdata = test.data, type = "response")
      
      try(cv.model <- FitGAM(data = train.data, reduce = FALSE, family.gam = gamfam, select = FALSE, link.fx = model$family$link, gam.formula = stats::formula(model), verbose = FALSE))
      
      if (exists("cv.model")) {
        cvpreds <- mgcv::predict.gam(object = cv.model, newdata = test.data, type = "response")
        
        if (strsplit(model$family$family, split = "[()]")[[1]][1] == "Negative Binomial") {
          cvtheta <- as.numeric(strsplit(cv.model$family[[1]], split = "[()]")[[1]][2])
          cvprobs <- 1 - stats::dnbinom(0, mu = cvpreds, size = cvtheta)
        } else {
          cvprobs <- 1 - stats::dpois(0, cvpreds)
        }
      } else {
        cvpreds <- rep(NA, times = nrow(test.data))
        cvprobs <- rep(NA, times = nrow(test.data))
        cv.model <- NA
      }
    }
    
    if (scale.preds) {
      scale.factors[i] <- mean(train.data[, species]) / mean(cvpreds)
      cvpreds <- cvpreds * scale.factors[i]
    }
    
    error.data$pred[start.vec[i]:end.vec[i]] <- preds
    error.data$prob[start.vec[i]:end.vec[i]] <- probs
    error.data$cvpred[start.vec[i]:end.vec[i]] <- cvpreds
    error.data$cvprob[start.vec[i]:end.vec[i]] <- cvprobs
    model.list[[i]] <- cv.model
    
    suppressWarnings(rm(cv.model, cvpreds))
    utils::setTxtProgressBar(pb, i)
  }
  
  if (scale.preds) {
    scale.factor <- mean(data[, species]) / mean(preds)
    error.data$pred <- error.data$pred * scale.factor
  }
  
  close(pb)
  error.data$error <- error.data$pred - error.data$abund
  error.data$cverror <- error.data$cvpred - error.data$abund
  
  if (scale.preds) {
    return(list(data = error.data, models = model.list, scale.factor = scale.factor, scale.factors = scale.factors))
  } else {
    return(list(data = error.data, models = model.list))
  }
}
