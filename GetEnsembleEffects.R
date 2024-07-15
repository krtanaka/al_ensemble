GetEnsembleEffects = function (effects.list, model.weights, vars = "all", scale = "log") 
{
  effects.list1 <- list()
  model.weights1 <- NULL
  list.index <- 1
  for (m in 1:length(effects.list)) {
    if (is.list(effects.list[[m]]) & model.weights[m] > 0) {
      effects.list1[[list.index]] <- effects.list[[m]]
      names(effects.list1)[list.index] <- names(effects.list)[m]
      model.weights1 <- c(model.weights1, model.weights[m])
      list.index <- list.index + 1
    }
  }
  var.table <- data.frame(type = NULL, dims = NULL, term = NULL)
  for (m in 1:length(effects.list1)) {
    for (n in 1:length(effects.list1[[m]])) {
      var.table0 <- data.frame(type = ifelse(nrow(effects.list1[[m]][[n]]) < 
                                               10, "factor", "smooth"), dims = ifelse("y" %in% 
                                                                                        names(effects.list1[[m]][[n]]), 2, 1), term = names(effects.list1[[m]])[n])
      var.table <- rbind(var.table, var.table0)
    }
  }
  var.table <- unique(var.table)
  if (vars[1] == "all") {
    var.table1 <- var.table
  }
  else {
    var.table1 <- var.table[var.table$term %in% vars]
  }
  vars2d <- var.table$term[var.table$dims == 2]
  vars1d <- var.table$term[var.table$dims == 1 & var.table$type == 
                             "smooth"]
  varsf <- var.table$term[var.table$type == "factor"]
  if ("lon*lat" %in% vars2d) {
    vars2d <- unique(c("lon*lat", vars2d))
  }
  out.list <- list()
  list.index <- 1
  if (length(vars2d) > 0) {
    for (v in 1:length(vars2d)) {
      temp <- effects.list1[[which(unlist(lapply(effects.list1, 
                                                 FUN = function(x) {
                                                   return(vars2d[v] %in% names(x))
                                                 })))[1]]]
      x.seq <- temp[[which(names(temp) == vars2d[v])]]$x
      y.seq <- temp[[which(names(temp) == vars2d[v])]]$y
      v2.preds <- as.data.frame(matrix(nrow = length(x.seq), 
                                       ncol = length(effects.list1)))
      colnames(v2.preds) <- names(model.weights1)
      for (m in 1:length(effects.list1)) {
        if (vars2d[v] %in% names(effects.list1[[m]])) {
          v2.preds[, m] <- exp(effects.list1[[m]][[vars2d[v]]]$effect) * 
            model.weights1[m]
        }
        else {
          v2.preds[, m] <- 0
        }
      }
      v2.dat <- data.frame(x = x.seq, y = y.seq, effect = apply(v2.preds, 
                                                                MARGIN = 1, FUN = sum))
      if (scale == "log") {
        v2.dat$effect[v2.dat$effect == 0] <- 0.01
        v2.dat$effect <- log(v2.dat$effect)
      }
      out.list[[list.index]] <- v2.dat
      names(out.list)[list.index] <- vars2d[v]
      list.index <- list.index + 1
    }
  }
  if (length(vars1d) > 0) {
    for (v in 1:length(vars1d)) {
      temp <- effects.list1[[which(unlist(lapply(effects.list1, 
                                                 FUN = function(x) {
                                                   return(vars1d[v] %in% names(x))
                                                 })))[1]]]
      x.seq <- temp[[which(names(temp) == vars1d[v])]]$x
      v.preds <- as.data.frame(matrix(nrow = length(x.seq), 
                                      ncol = length(model.weights1)))
      colnames(v.preds) <- names(model.weights1)
      m.vec <- NULL
      for (m in 1:length(effects.list1)) {
        if (vars1d[v] %in% names(effects.list1[[m]])) {
          v.preds[, m] <- exp(effects.list1[[m]][[vars1d[v]]]$effect) * 
            model.weights1[m]
          m.vec <- c(m.vec, m)
        }
        else {
          v.preds[, m] <- 0
        }
      }
      main.pred <- apply(X = v.preds, MARGIN = 1, FUN = sum)
      if (scale == "log") {
        main.pred[main.pred == 0] <- 0.01
        main.pred <- log(main.pred)
      }
      n.folds <- sum(unlist(lapply(strsplit(names(effects.list1[[m.vec[1]]][[vars1d[v]]]), 
                                            split = ""), FUN = function(x) {
                                              return(tolower(paste0(x[1:2], collapse = "")) == 
                                                       "cv")
                                            })))
      if (n.folds > 0) {
        cv.effects.dat <- as.data.frame(matrix(nrow = length(x.seq), 
                                               ncol = n.folds))
        for (f in 1:n.folds) {
          cv.preds <- as.data.frame(matrix(nrow = length(x.seq), 
                                           ncol = length(effects.list1)))
          for (m in 1:length(effects.list1)) {
            if (m %in% m.vec) {
              fold.match <- which(unlist(lapply(strsplit(names(effects.list1[[m]][[vars1d[v]]]), 
                                                         split = ""), FUN = function(x) {
                                                           return(tolower(paste0(x[1:2], collapse = "")) == 
                                                                    "cv")
                                                         })))
              cv.preds[, m] <- exp(effects.list1[[m]][[vars1d[v]]][, 
                                                                   fold.match[f]]) * model.weights1[m]
            }
            else {
              cv.preds[, m] <- 0
            }
          }
          cv.effects.dat[, f] <- apply(X = cv.preds, 
                                       MARGIN = 1, FUN = sum, na.rm = T)
          if (scale == "log") {
            cv.effects.dat[, f][cv.effects.dat[, f] == 
                                  0] <- 0.01
            cv.effects.dat[, f] <- log(cv.effects.dat[, 
                                                      f])
          }
        }
        names(cv.effects.dat) <- paste0("CV", 1:ncol(cv.effects.dat))
        uppers <- apply(X = cv.effects.dat, MARGIN = 1, 
                        FUN = "quantile", probs = 0.95, na.rm = T)
        lowers <- apply(X = cv.effects.dat, MARGIN = 1, 
                        FUN = "quantile", probs = 0.05, na.rm = T)
        out.list[[list.index]] <- data.frame(x = x.seq, 
                                             effect = main.pred, upper = uppers, lower = lowers, 
                                             cv.effects.dat)
        names(out.list)[list.index] <- vars1d[v]
      }
      else {
        out.list[[list.index]] <- data.frame(x.seq, effect = main.pred)
        names(out.list)[list.index] <- vars1d[v]
      }
      list.index <- list.index + 1
    }
  }
  if (length(varsf) > 0) {
    for (v in 1:length(varsf)) {
      temp <- effects.list1[[which(unlist(lapply(effects.list1, 
                                                 FUN = function(x) {
                                                   return(varsf[v] %in% names(x))
                                                 })))[1]]]
      f.seq <- factor(temp[[which(names(temp) == varsf[v])]]$x)
      f.preds <- as.data.frame(matrix(nrow = length(f.seq), 
                                      ncol = length(model.weights1)))
      m.vec <- NULL
      for (m in 1:length(model.weights1)) {
        if (varsf[v] %in% names(effects.list1[[m]])) {
          f.preds[, m] <- exp(effects.list1[[m]][[varsf[v]]]$effect) * 
            model.weights1[m]
          m.vec <- c(m.vec, m)
        }
        else {
          f.preds[, m] <- 0
        }
      }
      main.pred <- apply(X = f.preds, MARGIN = 1, FUN = sum)
      if (scale == "log") {
        main.pred[main.pred == 0] <- 0.01
        main.pred <- log(main.pred)
      }
      n.folds <- sum(unlist(lapply(strsplit(names(effects.list1[[m.vec[1]]][[varsf[v]]]), 
                                            split = ""), FUN = function(x) {
                                              return(tolower(paste0(x[1:2], collapse = "")) == 
                                                       "cv")
                                            })))
      if (n.folds > 0) {
        cv.effects.dat <- as.data.frame(matrix(nrow = length(f.seq), 
                                               ncol = n.folds))
        for (f in 1:n.folds) {
          cv.preds <- as.data.frame(matrix(nrow = length(f.seq), 
                                           ncol = length(effects.list1)))
          for (m in 1:length(effects.list1)) {
            if (m %in% m.vec) {
              fold.match <- which(unlist(lapply(strsplit(names(effects.list1[[m]][[varsf[v]]]), 
                                                         split = ""), FUN = function(x) {
                                                           return(tolower(paste0(x[1:2], collapse = "")) == 
                                                                    "cv")
                                                         })))
              cv.preds[, m] <- exp(effects.list1[[m]][[varsf[v]]][, 
                                                                  fold.match[f]]) * model.weights1[m]
            }
            else {
              cv.preds[, m] <- 0
            }
          }
          cv.effects.dat[, f] <- apply(X = cv.preds, 
                                       MARGIN = 1, FUN = sum, na.rm = T)
          if (scale == "log") {
            cv.effects.dat[, f][cv.effects.dat[, f] == 
                                  0] <- 0.01
            cv.effects.dat[, f] <- log(cv.effects.dat[, 
                                                      f])
          }
        }
        names(cv.effects.dat) <- paste0("CV", 1:ncol(cv.effects.dat))
        uppers <- apply(X = cv.effects.dat, MARGIN = 1, 
                        FUN = "quantile", probs = 0.95, na.rm = T)
        lowers <- apply(X = cv.effects.dat, MARGIN = 1, 
                        FUN = "quantile", probs = 0.05, na.rm = T)
        out.list[[list.index]] <- data.frame(x = f.seq, 
                                             effect = main.pred, upper = uppers, lower = lowers, 
                                             cv.effects.dat)
        names(out.list)[list.index] <- varsf[v]
      }
      else {
        out.list[[list.index]] <- data.frame(x = f.seq, 
                                             effect = main.pred)
        names(out.list)[list.index] <- varsf[v]
      }
      list.index <- list.index + 1
    }
  }
  return(out.list)
}