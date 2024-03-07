Effectsplot = function (effects.list, region = NA, crs = NA, nice.names = NULL, vars = "all") {
  
  vars <- vars[!is.na(vars)]
  
  if (length(vars) > 1 && all(vars != "all", na.rm = T) && all(is.character(vars))) {
    
    vars2 <- which(names(effects.list) %in% vars)
    
  } else {
    
    if (vars == "all") {
      vars2 <- 1:length(effects.list)
    } else {
      vars2 <- vars
    }
  }
  
  if (length(vars2) < length(vars)) {
    missing <- vars[which(vars %in% names(effects.list) == FALSE)]
    stop("Variables [", paste0(missing, collapse = ", "), "] not found in effects list")
  }
  
  effects.list2 <- list()
  
  for (i in 1:length(vars2)) {
    effects.list2[[i]] <- effects.list[[vars2[i]]]
    names(effects.list2)[i] <- names(effects.list)[vars2[i]]
  }
  
  if (is.null(nice.names)) {
    nice.names <- data.frame(var = unlist(strsplit(names(effects.list), "[*]")), 
                             name = unlist(strsplit(names(effects.list), "[*]")))
  }
  
  out.list <- list()
  
  for (j in 1:length(effects.list2)) {
    
    e.data <- effects.list2[[j]] %>% as.data.frame()
    
    if (nrow(e.data) == 1600) {
      
      effect.names <- strsplit(names(effects.list2)[j], split = "[*]")[[1]]
      xname <- nice.names$name[which(nice.names$var == effect.names[1])]
      yname <- nice.names$name[which(nice.names$var == effect.names[2])]
      
      if (xname %in% c("lon", "Longitude", "long")) {
        
        con.data <- grDevices::contourLines(x = sort(unique(e.data$x)), 
                                            y = sort(unique(e.data$y)), 
                                            z = matrix(nrow = 40, ncol = 40, data = e.data$effect))
        
        con.data2 <- data.frame(lon = con.data[[1]]$x, 
                                lat = con.data[[1]]$y, effect = con.data[[1]]$level, 
                                group = 1)
        label.spots <- data.frame(con.data2[round(nrow(con.data2) / 2), 1:2], 
                                  tag = con.data[[1]]$level)
        
        for (i in 2:length(con.data)) {
          c.dat <- data.frame(lon = con.data[[i]]$x, 
                              lat = con.data[[i]]$y, effect = con.data[[i]]$level, 
                              group = i)
          label.spots <- rbind(label.spots, 
                               data.frame(c.dat[round(nrow(c.dat) / 2), 1:2], 
                                          tag = con.data[[i]]$level))
          con.data2 <- rbind(con.data2, c.dat)
        }
        
        if (nrow(label.spots) > 10) {
          picks <- seq(from = 1, to = nrow(label.spots), 
                       by = floor(nrow(label.spots) / 10))
          label.spots <- label.spots[picks, ]
        }
        
        MAP <- akgfmaps::get_base_layers(select.region = tolower(region), 
                                         set.crs = "auto")
        if (is.na(crs)) {
          crs <- MAP$crs
        }
        
        e.ef <- sf::st_as_sf(x = con.data2, coords = c("lon", "lat"), crs = crs)
        e.ef2 <- sf::st_transform(e.ef, sf::st_crs(MAP$akland))
        spots.sf <- sf::st_as_sf(x = label.spots, coords = c("lon", "lat"), crs = crs)
        spots.sf2 <- sf::st_transform(spots.sf, sf::st_crs(MAP$akland))
        spot.data2 <- data.frame(sf::st_coordinates(spots.sf2), 
                                 label = spots.sf2$tag)
        e.data2 <- data.frame(sf::st_coordinates(e.ef2), 
                              e.ef2$effect, group = e.ef2$group)
        names(e.data2) <- c("lon", "lat", "effect", "group")
        ext.adjust.x <- c(0, 0)
        ext.adjust.y <- c(0, 0)
        
        
        if (tolower(region) == "goa") {
          ext.adjust.x <- c(4e+05, -2e+05)
          ext.adjust.y <- c(-490000, 6e+05)
          MAP$graticule$degree_label[c(1, 3, 5, 7, 9)] <- ""
          MAP$lon.breaks <- c(-170, -160, -150, -140, -130)
        }
        
        if (tolower(region) == "ai") {
          ext.adjust.x <- c(150000, -4e+05)
          ext.adjust.y <- c(-390000, 5e+05)
        }
        
        var.plot <- ggplot2::ggplot() + 
          ggplot2::geom_sf(data = MAP$akland, fill = "grey70") + 
          ggplot2::geom_sf(data = MAP$bathymetry, col = "grey60") + 
          ggplot2::geom_path(data = e.data2, ggplot2::aes(x = lon, y = lat, group = group), size = 1) + 
          ggplot2::geom_label(data = spot.data2, ggplot2::aes(x = X, y = Y, label = label), 
                              fill = grDevices::rgb(1, 1, 1, 0.9), label.size = NA, 
                              size = 4, label.padding = ggplot2::unit(0.1, "lines"), nudge_x = -1000) + 
          ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust.x, 
                            ylim = MAP$plot.boundary$y + ext.adjust.y) + 
          ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) + 
          ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) + 
          ggplot2::theme_bw() + 
          ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA), 
                         panel.background = ggplot2::element_rect(fill = NA, color = "black"), 
                         axis.title = ggplot2::element_blank(), 
                         axis.text = ggplot2::element_text(size = 12), 
                         plot.background = ggplot2::element_rect(fill = NA, color = NA))
        
      } else {
        
        con.data <- grDevices::contourLines(x = sort(unique(e.data[, 1])), 
                                            y = sort(unique(e.data[, 2])), 
                                            z = matrix(nrow = 40, ncol = 40, data = e.data$effect), nlevels = 10)
        
        con.data2 <- data.frame(x = con.data[[1]]$x, 
                                y = con.data[[1]]$y, effect = con.data[[1]]$level, 
                                group = 1)
        
        label.spots <- data.frame(con.data2[round(nrow(con.data2)/2), 
                                            1:2], tag = con.data[[1]]$level)
        
        for (i in 2:length(con.data)) {
          
          c.dat <- data.frame(x = con.data[[i]]$x, y = con.data[[i]]$y, effect = con.data[[i]]$level, group = i)
          label.spots <- rbind(label.spots, data.frame(c.dat[round(nrow(c.dat)/2), 1:2], tag = con.data[[i]]$level))
          con.data2 <- rbind(con.data2, c.dat)
          
        }
        
        names(e.data)[1:2] <- c("x", "y")
        
        var.plot <- ggplot2::ggplot() + ggplot2::geom_path(data = con.data2, 
                                                           ggplot2::aes(x = x, y = y, group = group)) + 
          ggplot2::geom_label(data = label.spots, ggplot2::aes(x = x, 
                                                               y = y, label = round(tag, 1)), fill = "white", 
                              label.size = NA) + ggplot2::xlab(xname) + 
          ggplot2::ylab(yname) + ggplot2::theme_bw() + 
          ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", 
                                                              fill = NA), panel.background = ggplot2::element_rect(fill = NA, 
                                                                                                                   color = "black"), axis.text = ggplot2::element_text(size = 12), 
                         axis.title = ggplot2::element_text(size = 12), 
                         plot.background = ggplot2::element_rect(fill = NA, 
                                                                 color = NA))
        if (xname == "Current Velocity East (m/s)") {
          var.plot <- var.plot + ggplot2::geom_vline(xintercept = 0, 
                                                     linetype = 3) + ggplot2::geom_hline(yintercept = 0, 
                                                                                         linetype = 3)
        }
        if (xname == "Current Velocity East SD (m/s)") {
          var.plot <- var.plot + ggplot2::geom_abline(intercept = 0, 
                                                      slope = 1, linetype = 3)
        }
      }
    }
    
    if (nrow(e.data) == 100) {
      
      xname <- nice.names$name[which(nice.names$var == names(effects.list2)[j])]
      span <- max(e.data$effect, na.rm = TRUE) - min(e.data$effect, na.rm = TRUE)
      
      upper.lim <- max(e.data$effect, na.rm = TRUE) + 3 * span
      lower.lim <- min(e.data$effect, na.rm = TRUE) - 3 * span
      
      y.lim <- c(ifelse(min(e.data$lower, na.rm = TRUE) < lower.lim, lower.lim, NA), 
                 ifelse(max(e.data$upper, na.rm = TRUE) > upper.lim, upper.lim, NA))
      
      var.plot <- ggplot2::ggplot() + 
        ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = effect)) + 
        ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = upper), linetype = 2) + 
        ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = lower), linetype = 2)
      
      if (any(!is.na(y.lim))) {
        var.plot <- var.plot + ggplot2::ylim(y.lim)
      }
      
      var.plot <- var.plot + 
        ggplot2::xlab(xname) + 
        ggplot2::ylab("Variable Effect") + 
        ggplot2::theme_bw() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA), 
                       panel.background = ggplot2::element_rect(fill = NA, color = "black"), 
                       axis.text = ggplot2::element_text(size = 12), 
                       axis.title = ggplot2::element_text(size = 12), 
                       plot.background = ggplot2::element_rect(fill = NA, color = NA))
      
    }
    
    if (nrow(e.data) < 10) {
      
      xname <- nice.names$name[which(nice.names$var == names(effects.list2)[j])]
      
      e.data$x <- as.numeric(as.character(e.data$x))
      
      var.plot <- ggplot2::ggplot() + 
        ggplot2::geom_segment(data = e.data, ggplot2::aes(y = effect, yend = effect, x = x - 0.35, xend = x + 0.35), size = 2) + 
        ggplot2::geom_segment(data = e.data, ggplot2::aes(y = lower, yend = lower, x = x - 0.35, xend = x + 0.35), size = 1, linetype = 2) + 
        ggplot2::geom_segment(data = e.data, ggplot2::aes(y = upper, yend = upper, x = x - 0.35, xend = x + 0.35), size = 1, linetype = 2) + 
        ggplot2::xlab(xname) + 
        ggplot2::ylab("Variable Effect") + 
        ggplot2::scale_x_continuous(breaks = (e.data$x)) + 
        ggplot2::theme_bw() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA), 
                       panel.background = ggplot2::element_rect(fill = NA, color = "black"), 
                       axis.text = ggplot2::element_text(size = 12), 
                       axis.title = ggplot2::element_text(size = 12), 
                       plot.background = ggplot2::element_rect(fill = NA, color = NA))
    }
    
    out.list[[j]] <- var.plot
    
    names(out.list)[j] <- names(effects.list2)[j]
    
  }
  
  return(out.list)
}