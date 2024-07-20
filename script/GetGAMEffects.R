GetGAMEffects = function (model, data, cv.model.list = NULL, vars = "all", scale = "log", 
                          scale.factor = 1) {
  
  if (is.null(cv.model.list) == F && all(is.na(cv.model.list))) {
    cv.model.list <- NULL
  }
  
  if (is.null(cv.model.list) == F) {
    list.index <- 1
    cv.model.list1 <- list()
    for (m in 1:length(cv.model.list)) {
      if (is.list(cv.model.list[[m]])) {
        cv.model.list1[[list.index]] <- cv.model.list[[m]]
        list.index <- list.index + 1
      }
    }
  }
  
  var.table <- AutodetectGAMTerms(model)
  
  if (model$family$family == "ziplss") {
    
    p.only <- var.table[[2]]$term[var.table[[2]]$term %in% var.table[[1]]$term == F]
    d.only <- var.table[[1]]$term[var.table[[1]]$term %in% var.table[[2]]$term == F]
    var.table <- unique(rbind(var.table[[1]], var.table[[2]]))
    
  } else {
    
    p.only <- NA
    
  }
  
  average.vars <- stats::na.omit(c(var.table$term[var.table$type %in% 
                                                    c("factor", "offset") == F], var.table$term2))
  
  fac.vars <- var.table$term[var.table$type == "factor"]
  off.var <- var.table$term[var.table$type == "offset"]
  
  if (tolower(vars) == "all") {
    do.rows <- which(var.table$type != "offset")
  }
  
  else {
    var.match <- unlist(strsplit(vars, split = "[*]"))
    do.rows <- vector(length = length(var.match))
    for (i in 1:length(var.match)) {
      do.rows[i] <- which(var.table$term == var.match[i] | 
                            var.table$term2 == var.match[i])
    }
    do.rows <- unique(do.rows)
  }
  
  average.dat <- apply(data.frame(data[, average.vars]), FUN = "mean", 
                       MARGIN = 2)
  
  fac.dat <- rep(0, times = length(fac.vars))
  off.dat <- apply(data.frame(data[, off.var]), FUN = "mean", 
                   MARGIN = 2)
  
  average.dat <- c(average.dat, fac.dat, off.dat)
  average.dat <- data.frame(t(average.dat))
  colnames(average.dat) <- c(average.vars, fac.vars, off.var)
  
  if (length(fac.vars) > 0) {
    
    for (f in 1:length(fac.vars)) {
      
      average.dat[, fac.vars[f]] <- factor(average.dat[, fac.vars[f]], levels = levels(as.factor(data[, fac.vars[f]])))
      
    }
  }
  
  effects.list <- list()
  
  list.index <- 1
  var.table1 <- var.table[do.rows, ]
  vars2d <- var.table1[var.table1$dims == 2, ]
  vars1d <- var.table1[var.table1$type == "smooth" & var.table1$dims == 1, ]
  varsf <- var.table1[var.table1$type == "factor", ]
  
  if (nrow(vars2d) > 0) {
    
    for (v in 1:nrow(vars2d)) {
      
      term2d <- c(vars2d$term[v], vars2d$term2[v])
      x.seq <- seq(from = min(data[, term2d[1]]), to = max(data[, 
                                                                term2d[1]]), length.out = 40)
      y.seq <- seq(from = min(data[, term2d[2]]), to = max(data[, 
                                                                term2d[2]]), length.out = 40)
      v2.data <- data.frame(x = rep(x.seq, times = 40), 
                            y = rep(y.seq, each = 40), average.dat[, names(average.dat) %in% 
                                                                     term2d == F])
      names(v2.data)[1:2] <- term2d
      v2.preds <- as.numeric(mgcv::predict.gam(model, type = "terms", 
                                               newdata = v2.data, terms = paste0("s(", term2d[1], 
                                                                                 ",", term2d[2], ")"))[, 1])
      if (model$family$family == "ziplss") {
        
        v2.probs <- as.numeric(mgcv::predict.gam(model, type = "terms", newdata = v2.data, terms = paste0("s.1(", term2d[1], ",", term2d[2], ")"))[, 1])
        v2.preds <- exp(v2.preds) * (1 - exp(-exp(v2.probs)))/(1 - stats::dpois(0, exp(v2.preds)))
        v2.preds <- log(v2.preds)
        
      }
      
      if (scale == "abund") {
        v2.preds <- exp(v2.preds) * scale.factor
      } else {
        v2.preds <- v2.preds + log(scale.factor)
      }
      
      grDevices::png("trashme.png")
      x <- plot(model, scale = 0, se = F, pages = 1)
      grDevices::dev.off()
      file.remove("trashme.png")
      xvec <- vector(length = length(x))
      for (n in 1:length(x)) {
        xvec[n] <- x[[n]]$xlab == term2d[1]
      }
      v2nas <- which(is.na(x[[which(xvec)[1]]]$fit))
      v2.preds[v2nas] <- NA
      v2.preds[is.infinite(v2.preds)] <- NA
      effects.list[[list.index]] <- data.frame(x.seq, y.seq, 
                                               effect = v2.preds)
      names(effects.list[[list.index]])[1:2] <- c("x", 
                                                  "y")
      names(effects.list)[list.index] <- paste(term2d, 
                                               collapse = "*")
      list.index <- list.index + 1
    }
  }
  
  if (nrow(vars1d) > 0) {
    for (v in 1:nrow(vars1d)) {
      term1d <- vars1d$term[v]
      v.data <- data.frame(seq(from = min(data[, term1d]), 
                               to = max(data[term1d]), length.out = 100), average.dat[names(average.dat) != 
                                                                                        term1d])
      names(v.data)[1] <- term1d
      v.name <- paste0("s(", term1d, ")")
      main.pred <- as.numeric(mgcv::predict.gam(model, 
                                                newdata = v.data, type = "terms", terms = v.name))
      se.pred <- as.numeric(mgcv::predict.gam(model, newdata = v.data, 
                                              type = "terms", terms = v.name, se.fit = T)[[2]])
      if (term1d %in% p.only) {
        main.pred <- main.pred[1:nrow(v.data)]
        se.pred <- se.pred[1:nrow(v.data)]
      }
      if (model$family$family == "ziplss") {
        main.pred[-10 > main.pred] <- log(0.01)
        main.prob <- as.numeric(mgcv::predict.gam(model, 
                                                  type = "terms", newdata = v.data, terms = paste0("s.1(", 
                                                                                                   term1d, ")"))[, 1])
        main.pred <- exp(main.pred) * (1 - exp(-exp(main.prob)))/(1 - 
                                                                    stats::dpois(0, exp(main.pred)))
        main.pred <- log(main.pred)
      }
      if (scale == "abund") {
        main.pred <- exp(main.pred) * scale.factor
      }
      else {
        main.pred <- main.pred + log(scale.factor)
      }
      main.pred[is.infinite(main.pred)] <- NA
      if (is.null(cv.model.list) == F) {
        effect.dat <- matrix(nrow = 100, ncol = length(cv.model.list))
        for (f in 1:length(cv.model.list)) {
          if (is.list(cv.model.list[[f]])) {
            cv.pred <- as.numeric(mgcv::predict.gam(cv.model.list[[f]], 
                                                    newdata = v.data, type = "terms", terms = v.name))
            if (term1d %in% p.only) {
              cv.pred <- cv.pred[1:nrow(v.data)]
            }
            if (model$family$family == "ziplss") {
              cv.pred[-10 > cv.pred] <- log(0.01)
              cv.prob <- as.numeric(mgcv::predict.gam(cv.model.list[[f]], 
                                                      type = "terms", newdata = v.data, terms = paste0("s.1(", 
                                                                                                       term1d, ")"))[, 1])
              cv.pred <- exp(cv.pred) * (1 - exp(-exp(cv.prob)))/(1 - 
                                                                    stats::dpois(0, exp(cv.pred)))
              cv.pred <- log(cv.pred)
            }
            if (scale == "abund") {
              cv.pred <- exp(cv.pred) * scale.factor
            }
            else {
              cv.pred <- cv.pred + log(scale.factor)
            }
          }
          else {
            cv.pred <- NA
          }
          effect.dat[, f] <- cv.pred
        }
        colnames(effect.dat) <- paste0("CV", 1:length(cv.model.list))
        effect.dat[is.infinite(effect.dat)] <- NA
        uppers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", 
                        probs = 0.95, na.rm = T)
        lowers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", 
                        probs = 0.05, na.rm = T)
        variance <- apply(X = effect.dat, MARGIN = 1, 
                          FUN = var, na.rm = T)
        out.dat <- data.frame(x = v.data[, 1], effect = main.pred, 
                              var = se.pred^2, cvvar = variance, upper = uppers, 
                              lower = lowers, effect.dat)
      }
      else {
        out.dat <- data.frame(x = v.data[, 1], effect = main.pred, 
                              var = se.pred^2)
      }
      effects.list[[list.index]] <- out.dat
      names(effects.list)[list.index] <- term1d
      list.index <- list.index + 1
    }
  }
  if (nrow(varsf) > 0) {
    for (v in 1:nrow(varsf)) {
      termf <- varsf$term[v]
      if (termf %in% names(model$model)) {
        fac.name <- termf
      }
      else {
        fac.name <- paste0("as.factor(", termf, ")")
      }
      f.data <- data.frame(unique(as.character(data[, termf])), 
                           average.dat[names(average.dat) != termf])
      names(f.data)[1] <- termf
      main.pred <- as.numeric(mgcv::predict.gam(model, 
                                                newdata = f.data, type = "terms", terms = fac.name))
      se.pred <- as.numeric(mgcv::predict.gam(model, newdata = f.data, 
                                              type = "terms", terms = fac.name, se.fit = T)[[2]])
      if (termf %in% p.only) {
        main.pred <- main.pred[1:nrow(f.data)]
        se.pred <- se.pred[1:nrow(f.data)]
      }
      if (model$family$family == "ziplss") {
        main.pred[-10 > main.pred] <- log(0.01)
        main.prob <- as.numeric(mgcv::predict.gam(model, 
                                                  type = "terms", newdata = f.data, terms = paste0("as.factor(", 
                                                                                                   termf, ").1"))[, 1])
        main.pred <- exp(main.pred) * (1 - exp(-exp(main.prob)))/(1 - 
                                                                    stats::dpois(0, exp(main.pred)))
        main.pred <- log(main.pred)
      }
      if (scale == "abund") {
        main.pred <- exp(main.pred) * scale.factor
      }
      else {
        main.pred <- main.pred + log(scale.factor)
      }
      if (is.null(cv.model.list) == F) {
        effect.dat <- matrix(nrow = length(main.pred), 
                             ncol = length(cv.model.list))
        for (f in 1:length(cv.model.list)) {
          if (is.list(cv.model.list[[f]])) {
            cv.pred <- as.numeric(mgcv::predict.gam(cv.model.list[[f]], 
                                                    newdata = f.data, type = "terms", terms = fac.name))
            if (termf %in% p.only) {
              cv.pred <- cv.pred[1:nrow(f.data)]
            }
            if (model$family$family == "ziplss") {
              cv.pred[-10 > cv.pred] <- log(0.01)
              cv.prob <- as.numeric(mgcv::predict.gam(cv.model.list[[f]], 
                                                      type = "terms", newdata = f.data, terms = paste0("as.factor(", 
                                                                                                       termf, ").1"))[, 1])
              cv.pred <- exp(cv.pred) * (1 - exp(-exp(cv.prob)))/(1 - 
                                                                    stats::dpois(0, exp(cv.pred)))
              cv.pred <- log(cv.pred)
            }
          }
          else {
            cv.pred <- NA
          }
          if (scale == "abund") {
            cv.pred <- exp(cv.pred) * scale.factor
          }
          else {
            cv.pred <- cv.pred + log(scale.factor)
          }
          effect.dat[, f] <- cv.pred
        }
        colnames(effect.dat) <- paste0("CV", 1:length(cv.model.list))
        uppers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", 
                        probs = 0.95, na.rm = T)
        lowers <- apply(X = effect.dat, MARGIN = 1, FUN = "quantile", 
                        probs = 0.05, na.rm = T)
        variance <- apply(X = effect.dat, MARGIN = 1, 
                          FUN = var, na.rm = T)
        out.dat <- data.frame(x = f.data[, 1], effect = main.pred, 
                              var = se.pred^2, cvvar = variance, upper = uppers, 
                              lower = lowers, effect.dat)
      }
      else {
        out.dat <- data.frame(x = f.data[, 1], effect = main.pred, 
                              var = se.pred^2)
      }
      effects.list[[list.index]] <- out.dat
      names(effects.list)[list.index] <- termf
      list.index <- list.index + 1
    }
  }
  return(effects.list)
}