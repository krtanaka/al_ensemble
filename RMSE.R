RMSE = function (pred, obs) {
  keep <- which(is.na(pred) == F & is.na(obs) == F)
  return(sqrt(sum((pred[keep] - obs[keep])^2)/length(pred[keep])))
}