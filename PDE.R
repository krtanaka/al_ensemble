PDE = function (obs, pred) 
{
    term1 <- obs * log(obs/mean(obs))
    term1[is.nan(term1)] <- 0
    term2 <- obs - mean(obs)
    nulldev <- 2 * sum(term1 - term2)
    pred[pred < 1e-05] <- 1e-05
    term1 <- obs * log(obs/pred)
    term1[is.nan(term1)] <- 0
    term2 <- obs - pred
    pdev <- 2 * sum(term1 - term2)
    return(1 - (pdev/nulldev))
}