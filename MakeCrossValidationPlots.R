MakeCrossValidationPlots = function (error.data, method = "pearson", make.hist = TRUE) {
  
  keepers <- unique(which(is.infinite(error.data$pred) == FALSE & 
                            is.infinite(error.data$cvpred) == FALSE))
  
  if (length(keepers) < nrow(error.data)) {
    warning(paste0(nrow(error.data) - length(keepers), 
                   " Infinite values detected and removed; estimates of error may be too low"))
  }
  
  error.data2 <- error.data
  
  if (method == "spearman") {
    error.data2$abund <- rank(error.data2$abund)
    error.data2$pred <- rank(error.data2$pred)
    error.data2$cvpred <- rank(error.data2$cvpred)
  }
  
  main.regr <- stats::lm(error.data2$pred[keepers] ~ error.data2$abund[keepers])
  main.r2 <- summary(main.regr)$r.squared
  main.rmse <- sqrt(sum((stats::na.omit(error.data$abund[keepers] - 
                                          error.data$pred[keepers]))^2) / 
                      nrow(stats::na.omit(error.data[keepers, ])))
  
  cv.regr <- stats::lm(error.data2$cvpred[keepers] ~ error.data2$abund[keepers])
  cv.r2 <- summary(cv.regr)$r.squared
  cv.rmse <- sqrt(sum((stats::na.omit(error.data$abund[keepers] - 
                                        error.data$cvpred[keepers]))^2) / 
                    nrow(stats::na.omit(error.data[keepers, ])))
  
  print(paste("Full model", method, "Rsq =", round(main.r2, 2)))
  print(paste("CV", method, "Rsq =", round(cv.r2, 2)))
  print(paste("Full model RMSE =", round(main.rmse, 3)))
  print(paste("CV RMSE =", round(cv.rmse, 3)))
  
  old.par <- graphics::par()[c("mfcol", "family", "mar", "xaxs", "yaxs")]
  
  graphics::par(mfcol = c(ifelse(make.hist == TRUE, 3, 2), 2), 
                family = "sans", mar = c(4, 4, 3, 1))
  
  stats::qqnorm((error.data$pred[keepers] - error.data$abund[keepers]), main = "Model Predictions", bty = "l")
  stats::qqline((error.data$pred[keepers] - error.data$abund[keepers]))
  
  if (make.hist == TRUE) {
    graphics::hist((error.data$pred[keepers] - error.data$abund[keepers]), 
                   xlab = "Residuals", main = "")
  }
  
  pred.max <- ifelse(method == "pearson", 
                     stats::quantile(error.data2$pred[keepers], probs = 0.99, na.rm = TRUE), 
                     nrow(error.data2))
  
  abund.max <- ifelse(method == "pearson", 
                      stats::quantile(error.data2$pred[keepers], probs = 0.99, na.rm = TRUE), 
                      nrow(error.data2))
  
  plot.max <- max(pred.max, abund.max) * 1.1
  
  plot(y = error.data2$pred[keepers], x = error.data2$abund[keepers], 
       ylim = c(0, plot.max), xlim = c(0, plot.max), 
       ylab = ifelse(method == "pearson", "Predicted", "Predicted Ranks"), 
       xlab = ifelse(method == "pearson", "Observed", "Observed Ranks"), 
       main = "", pch = 20, bty = "l")
  graphics::abline(coef = c(0, 1), lty = 2)
  graphics::abline(main.regr, col = 2)
  graphics::text(1, plot.max * 0.9, paste(method, "R-squared = ", signif(main.r2, 2)), pos = 4)
  
  pred.max <- ifelse(method == "pearson", 
                     stats::quantile(error.data2$cvpred[keepers], probs = 0.99, na.rm = TRUE), 
                     nrow(error.data2))
  abund.max <- ifelse(method == "pearson", 
                      stats::quantile(error.data2$cvpred[keepers], probs = 0.99, na.rm = TRUE), 
                      nrow(error.data2))
  plot.max <- max(pred.max, abund.max) * 1.1
  
  stats::qqnorm((error.data$cvpred[keepers] - error.data$abund[keepers]), main = "Test data", bty = "l")
  stats::qqline((error.data$cvpred[keepers] - error.data$abund[keepers]))
  
  if (make.hist == TRUE) {
    graphics::hist((error.data$cvpred[keepers] - error.data$abund[keepers]), 
                   xlab = "Residuals", main = "")
  }
  
  plot(y = error.data2$cvpred[keepers], x = error.data2$abund[keepers], 
       ylim = c(0, plot.max), xlim = c(0, plot.max), 
       ylab = ifelse(method == "pearson", "Predicted", "Predicted Ranks"), 
       xlab = ifelse(method == "pearson", "Observed", "Observed Ranks"), 
       main = "", pch = 20, bty = "l")
  graphics::abline(coef = c(0, 1), lty = 2)
  graphics::abline(cv.regr, col = 2)
  graphics::text(1, plot.max * 0.9, paste(method, "R-squared = ", signif(cv.r2, 2)), pos = 4)
  
  suppressWarnings(graphics::par(old.par))
}
