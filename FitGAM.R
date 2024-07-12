FitGAM = function (data, gam.formula = NULL, reduce = F, select = F, verbose = F, 
                   family.gam = "poisson", theta = NA, link.fx = NA) 
{
  gam.form <- stats::as.formula(gam.formula)
  species <- as.character(gam.form)[[2]]
  off.0 <- trimws(strsplit(as.character(gam.form)[[3]], split = "[+]")[[1]])
  off.1 <- unlist(lapply(strsplit(off.0, split = "[()]"), FUN = function(x) {
    return(ifelse(x[1] == "offset", x[2], NA))
  }))
  offset.gam <- off.1[is.na(off.1) == F]
  if (length(offset.gam) == 0) {
    offset.gam <- NULL
  }
  if (family.gam %in% c("gaussian", "quasigaussian") & is.na(link.fx)) {
    link.fx <- "identity"
  }
  if (family.gam %in% c("poisson", "quasipoisson", "nb") & 
      is.na(link.fx)) {
    link.fx <- "log"
  }
  if (family.gam == "binomial") {
    if (is.na(link.fx)) {
      link.fx <- "logit"
    }
    data[, species] <- as.integer(data[, species] > 0)
  }
  if (is.na(theta) == F) {
    gam.fam <- eval(paste0(family.gam, "(theta=", theta, 
                           ",link=", link.fx, ")"))
  }
  else {
    gam.fam <- eval(paste0(family.gam, "(link=", link.fx, 
                           ")"))
  }
  xvar.gam <- trimws(unlist(strsplit(as.character(gam.form)[[3]], 
                                     split = "[+]")))
  for (v in 1:length(xvar.gam)) {
    if (strsplit(xvar.gam[v], split = "[(]")[[1]][1] == "offset") {
      xvar.gam <- xvar.gam[-v]
      break
    }
  }
  if (is.null(offset.gam) == F) {
    offset.form <- paste0("offset(", offset.gam, ")")
  }
  else {
    offset.form <- NULL
  }
  for (i in 1:length(xvar.gam)) {
    phrase <- strsplit(xvar.gam, split = " ")[[i]]
    spot <- which(phrase == "bs") + 2
    if (length(spot) > 0) {
      phrase2 <- strsplit(phrase[spot], "")[[1]]
      lets <- which(phrase2 %in% letters)
      newphrase <- paste0("'", paste(phrase2[lets], collapse = ""), 
                          "',")
      phrase[spot] <- newphrase
      xvar.gam[i] <- paste(phrase, collapse = "")
    }
  }
  if (select) {
    out.gam <- mgcv::gam(gam.form, family = gam.fam, data = data, 
                         select = T)
    if (verbose) {
      print(summary(out.gam))
    }
    worst.var <- which.min(summary(out.gam)$edf)
    worst.name <- names(summary(out.gam)$s.table[, 1])[worst.var]
    target.edf <- length(unlist(strsplit(worst.name, split = "[(,)]"))) - 
      1
    if (summary(out.gam)$edf[worst.var] < 1) {
      bad.vars <- worst.var
    }
    if (summary(out.gam)$edf[worst.var] >= 1) {
      bad.vars <- NA
    }
    while (is.na(bad.vars) == F) {
      xvar.gam <- xvar.gam[-bad.vars]
      gam.form <- stats::as.formula(paste0(paste(species, 
                                                 "~", paste(c(xvar.gam, offset.form), collapse = "+"))))
      out.gam <- mgcv::gam(gam.form, family = gam.fam, 
                           data = data, select = T)
      worst.var <- which.min(summary(out.gam)$edf)
      if (summary(out.gam)$edf[worst.var] < 1) {
        bad.vars <- worst.var
      }
      if (summary(out.gam)$edf[worst.var] >= 1) {
        bad.vars <- NA
      }
      if (verbose) {
        print(summary(out.gam))
      }
    }
  }
  else {
    out.gam <- mgcv::gam(gam.form, family = gam.fam, data = data, 
                         select = F)
  }
  if (reduce) {
    while (reduce) {
      gcv_gam <- out.gam$gcv.ubre
      pvals <- summary(out.gam)$s.pv
      pvals <- c(pvals, summary(out.gam)$p.pv[-1])
      least_sig <- which.max(pvals)
      xvar.gam1 <- xvar.gam[-least_sig]
      gam.form1 <- stats::as.formula(paste0(paste(species, 
                                                  "~", paste(c(xvar.gam1, offset.form), collapse = "+"))))
      gam.form2.1 <- stats::as.formula(paste0(as.character(gam.form1)[1], 
                                              as.character(gam.form1)[3]))
      out.gam1 <- mgcv::gam(gam.form1, family = gam.fam, 
                            data = data)
      gcv_gam1 <- out.gam1$gcv.ubre
      if (gcv_gam > gcv_gam1) {
        gam.form <- gam.form1
        gcv_gam <- gcv_gam1
        xvar.gam <- xvar.gam1
        out.gam <- out.gam1
      }
      if (verbose == T) {
        print(summary(out.gam))
      }
      if (gcv_gam <= gcv_gam1) {
        reduce <- F
      }
    }
  }
  return(out.gam)
}