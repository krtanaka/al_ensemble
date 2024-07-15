AutodetectGAMTerms = function (model, hgam = "all") {
  
  n.formulas <- 1
  
  if (model$family$family == "ziplss" & hgam %in% c("b", "both", "all")) {
    n.formulas <- 2
  }
  
  for (f in 1:n.formulas) {
    
    form.index <- 3
    
    if (model$family$family == "ziplss") {
      if (hgam %in% c("b", "both", "all")) {
        form1 <- stats::formula(model)[[f]]
      }
      if (hgam %in% c("d", "dens", "density")) {
        form1 <- stats::formula(model)[[1]]
      }
      if (hgam %in% c("p", "prob", "probability")) {
        form1 <- stats::formula(model)[[2]]
      }
    } else {
      form1 <- stats::formula(model)
    }
    
    terms <- trimws(strsplit(as.character(form1[[length(form1)]]), split = "[+]")[[2]])
    type.dat <- data.frame(
      type = rep(NA, length(terms)), 
      dims = 1, 
      term = NA, 
      term2 = NA, 
      bs = NA, 
      k = NA, 
      m = NA, 
      m2 = NA
    )
    
    for (t in 1:length(terms)) {
      x <- terms[t]
      x2 <- strsplit(x, split = "[(=)]")[[1]]
      
      if (length(x2) == 1) {
        type.dat$type[t] <- "linear"
        type.dat$term[t] <- x2
      } else {
        if (x2[1] %in% c("s", "te")) {
          dims <- length(strsplit(x2[2], split = ", ")[[1]]) - 1
          for (n in 1:dims) {
            type.dat[t, 2 + n] <- strsplit(x2[2], split = ", ")[[1]][n]
          }
          type.dat$type[t] <- "smooth"
          type.dat$dims[t] <- dims
          
          formula.options <- strsplit(x, split = ",")[[1]]
          
          bs.spot <- which(unlist(lapply(strsplit(formula.options, "="), FUN = function(x) {
            return(trimws(x[1]))
          })) == "bs")
          
          if (length(bs.spot) > 0) {
            type.dat$bs[t] <- strsplit(formula.options[bs.spot], split = "\"")[[1]][2]
          }
          
          k.spot <- which(unlist(lapply(strsplit(formula.options, "="), FUN = function(x) {
            return(trimws(x[1]))
          })) == "k")
          
          if (length(k.spot) > 0) {
            type.dat$k[t] <- trimws(strsplit(strsplit(formula.options[k.spot], split = "=")[[1]][2], split = "[)]")[[1]])
          }
          
          m.spot <- which(unlist(lapply(strsplit(formula.options, "="), FUN = function(x) {
            return(trimws(x[1]))
          })) == "m")
          
          if (length(m.spot) > 0) {
            if ("c" %in% strsplit(formula.options[m.spot], split = "")[[1]]) {
              m.spot <- c(m.spot, m.spot + 1)
            }
            for (n in 1:length(m.spot)) {
              m1 <- trimws(formula.options[m.spot][n])
              if (n == 1) {
                m2 <- trimws(strsplit(m1, split = "=")[[1]][2])
              } else {
                m2 <- m1
              }
              m3 <- trimws(strsplit(m2, split = "[()]")[[1]])
              type.dat[t, 6 + n] <- suppressWarnings(stats::na.omit(as.numeric(m3))[1])
            }
          }
        }
        
        if (x2[1] == "as.factor") {
          type.dat$type[t] <- "factor"
          type.dat$term[t] <- x2[2]
        }
      }
    }
    
    terms2 <- unlist(strsplit(x = names(model$model), split = "[()]"))
    off.term <- which(terms2 == "offset") + 1
    
    if (length(off.term) > 0) {
      type.dat <- rbind(type.dat, data.frame(
        type = "offset", 
        dims = 1, 
        term = terms2[off.term], 
        term2 = NA, 
        bs = NA, 
        k = NA, 
        m = NA, 
        m2 = NA
      ))
    }
    
    if (n.formulas > 1) {
      if (f == 1) {
        out.dat <- list(type.dat)
      } else {
        out.dat[[f]] <- type.dat
      }
    } else {
      out.dat <- type.dat
    }
  }
  
  return(out.dat)
}