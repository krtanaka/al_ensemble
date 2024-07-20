MakeAKGFDotplot = function (presence, absence = NA, highdensity = NA, dataCRS = NA, 
                            region, survey.area = "default", ext.adjust = "default", 
                            legend.pos = "default", title.name = NA, title.count = T, 
                            title.pos = "default", abs.name = "absent", pres.name = "present", 
                            hd.name = "top 10%", abs.col = "steelblue4", pres.col = "orange", 
                            hd.col = "red", abs.shape = 16, pres.shape = 1, hd.shape = 16, 
                            abs.size = 0.1, pres.size = 1, hd.size = 1) {
  
  presence = species.data[species.data[, s] > 0, ]
  absence = species.data[species.data[, s] == 0, ]
  highdensity = species.data[species.data[, s] > hd, ]
  region = tolower(region)
  dataCRS = raster.stack@crs
  title.name = figure.name
  
  region <- tolower(region)
  
  if (region %in% c("ebs", "bs.all", "sebs", "bs.south", "ecs", 
                    "ebs.ecs", "ai", "ai.west", "ai.central", "ai.east", 
                    "goa", "goa.west", "goa.east")) {
    MAP <- akgfmaps::get_base_layers(select.region = region, 
                                     set.crs = "auto", use.survey.bathymetry = TRUE)
  } else {
    stop("region not recognized")
  }
  
  if (ext.adjust == "default") {
    if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("ebs", "bs.all", "bs.south", "sebs", 
                      "ecs", "ebs.ecs")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("goa", "goa.west", "goa.east")) {
      ext.adjust <- c(2e+05, -1e+05, 0, 0)
    }
  }
  
  if (length(legend.pos) < 2) {
    
    if (legend.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", 
                        "ai.east")) {
        legend.pos <- c(0.12, 0.2)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", 
                        "ecs", "ebs.ecs")) {
        legend.pos <- c(0.12, 0.18)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        legend.pos <- c(0.08, 0.41)
      }
    } else {
      legend.pos <- NA
    }
  }
  
  if (length(title.pos) < 2) {
    
    if (title.pos == "default") {
      
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        title.pos <- c(-9e+05, 490000)
      }
      
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        title.pos <- c(-550000, 8e+05)
      }
      
      if (region %in% c("goa", "goa.west", "goa.east")) {
        title.pos <- c(-1300000, 630000)
      }
    } else {
      title.pos <- NA
    }
  }
  
  if (is.character(survey.area) && survey.area[1] == "default") {
    
    survey.sf <- MAP$survey.area
    
  }else {
    
    if (class(survey.area)[1] == "sf") {
      survey.sf <- sf::st_transform(survey.area, sf::st_crs(MAP$akland))
    }
    
    if (class(survey.area)[1] == "RasterLayer") {
      survey.sf1 <- stars::st_as_stars(is.na(survey.area))
      survey.sf2 <- sf::st_as_sf(survey.sf1, merge = TRUE)
      survey.sf3 <- sf::st_cast(survey.sf2, "POLYGON")
      survey.sf <- sf::st_transform(survey.sf3, sf::st_crs(MAP$akland))[1:(nrow(survey.sf3) - 1), ]
    }
  }
  
  leg.col <- NULL
  leg.shape <- NULL
  leg.size <- NULL
  leg.name <- NULL
  
  if (is.data.frame(absence)) {
    
    if (is.na(dataCRS)) {
      
      abs.sf <- sf::st_as_sf(x = absence, coords = c("lon", "lat"), crs = sf::st_crs(MAP$akland))
      
    } else {
      
      abs.sf0 <- sf::st_as_sf(x = absence, coords = c("lon", "lat"), crs = dataCRS)
      abs.sf <- sf::st_transform(abs.sf0, sf::st_crs(MAP$akland))
      
    }
    
    leg.name <- abs.name
    leg.col <- abs.col
    leg.shape <- abs.shape
    leg.size <- abs.size
    abs.fac <- 1
    
  } else {
    
    abs.fac <- 0
    
  }
  
  if (is.na(dataCRS)) {
    
    pres.sf <- sf::st_as_sf(x = presence, coords = c("lon", "lat"), crs = sf::st_crs(MAP$akland))
    
  } else {
    
    pres.sf0 <- sf::st_as_sf(x = presence, coords = c("lon", "lat"), crs = dataCRS)
    pres.sf <- sf::st_transform(pres.sf0, sf::st_crs(MAP$akland))
    
  }
  
  leg.name <- c(leg.name, pres.name)
  leg.col <- c(leg.col, pres.col)
  leg.shape <- c(leg.shape, pres.shape)
  leg.size <- c(leg.size, pres.size)
  pres.fac <- abs.fac + 1
  
  if (is.data.frame(highdensity)) {
    
    if (is.na(dataCRS)) {
      high.sf <- sf::st_as_sf(x = highdensity, coords = c("lon", "lat"), crs = sf::st_crs(MAP$akland))
      
    } else {
      high.sf0 <- sf::st_as_sf(x = highdensity, coords = c("lon", "lat"), crs = dataCRS)
      high.sf <- sf::st_transform(high.sf0, sf::st_crs(MAP$akland))
    }
    
    leg.name <- c(leg.name, hd.name)
    leg.col <- c(leg.col, hd.col)
    leg.shape <- c(leg.shape, hd.shape)
    leg.size <- c(leg.size, hd.size)
    hd.fac <- pres.fac + 1
  }
  
  dotplot <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = survey.sf, fill = "grey95") + 
    # ggplot2::geom_sf(data = MAP$akland, fill = "grey40") + 
    # ggplot2::geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) + 
    ggplot2::geom_sf(data = MAP$bathymetry, color = "grey60")
  
  if (is.data.frame(absence)) {
    
    dotplot <-
      # dotplot +
      ggplot2::ggplot() + 
      ggplot2::geom_sf(data = abs.sf, 
                       alpha = 0.25, size = abs.size, shape = abs.shape, 
                       ggplot2::aes(color = factor(abs.fac)))
  }
  
  dotplot <- dotplot + ggplot2::geom_sf(data = pres.sf, size = pres.size, 
                                        ggplot2::aes(color = factor(pres.fac)), 
                                        shape = pres.shape, 
                                        stroke = 0.8)
  if (is.data.frame(highdensity)) {
    dotplot <- dotplot + 
      ggplot2::geom_sf(data = high.sf, size = hd.size, shape = hd.shape, 
                       ggplot2::aes(color = factor(hd.fac)))
  }
  
  dotplot <- dotplot + 
    # ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust[1:2], 
    # ylim = MAP$plot.boundary$y + ext.adjust[3:4]) + 
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) + 
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA), 
                   panel.background = ggplot2::element_rect(fill = NA, color = "black"),
                   legend.key = ggplot2::element_rect(fill = NA, color = "grey30"), 
                   legend.position = legend.pos, 
                   axis.title = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_text(size = 12), 
                   legend.text = ggplot2::element_text(size = 12), 
                   legend.title = ggplot2::element_text(size = 12), 
                   plot.background = ggplot2::element_rect(fill = NA, color = NA))
  
  if (is.na(title.name) == F && length(title.pos) == 2) {
    
    if (title.count) {
      
      title.name <- paste0(title.name, "\nN = ", format(nrow(pres.sf), big.mark = ","))
      
    }
    
    dotplot <- dotplot + 
      ggplot2::geom_label(data = data.frame(x = title.pos[1], y = title.pos[2], label = title.name), 
                          ggplot2::aes(x = x, y = y, label = label, hjust = 0, vjust = 1), size = 5)
  }
  
  if (length(legend.pos) == 2) {
    
    dotplot <- dotplot + 
      ggplot2::scale_color_manual(name = NULL, values = leg.col, labels = leg.name) + 
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = leg.shape, size = leg.size)))
    
  }
  
  return(dotplot)
}