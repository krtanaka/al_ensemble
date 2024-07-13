MakeAKGFDensityplot = function (
    region, 
    density.map,
    buffer = 1, 
    survey.area = "default", 
    ext.adjust = "default", 
    legend.pos = "default", 
    legend.title = NA, 
    col.palette = "plasma",
    col.palette.limits = c(0, 1),
    title.name = NA, 
    title.pos = "default", 
    barheight = 5) {
  
  region <- tolower(region)
  
  if (region %in% c("ebs", "bs.all", "sebs", "bs.south", "ecs", 
                    "ebs.ecs", "ai", "ai.west", "ai.central", "ai.east", 
                    "goa", "goa.west", "goa.east")) {
    
    MAP <- akgfmaps::get_base_layers(select.region = region, set.crs = "auto")
    
  } else {
    stop("region not recognized")
  }
  
  if (ext.adjust == "default") {
    if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("goa", "goa.west", "goa.east")) {
      ext.adjust <- c(200000, -100000, 0, 0)
    }
  }
  
  if (length(legend.pos) < 2) {
    if (legend.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        legend.pos <- c(0.08, 0.28)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        legend.pos <- c(0.12, 0.18)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        legend.pos <- c(0.08, 0.51)
      }
    } else {
      legend.pos <- NA
    }
  }
  
  if (length(title.pos) < 2) {
    if (title.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        title.pos <- c(-900000, 490000)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        title.pos <- c(-550000, 800000)
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
  } else {
    if (class(survey.area)[1] == "sf") {
      survey.sf <- sf::st_transform(survey.area, sf::st_crs(MAP$akland))[1:(nrow(survey.area) - 1), ]
    }
    if (class(survey.area)[1] == "RasterLayer") {
      survey.sf1 <- stars::st_as_stars(is.na(survey.area))
      survey.sf2 <- sf::st_as_sf(survey.sf1, merge = TRUE)
      survey.sf3 <- sf::st_cast(survey.sf2, "POLYGON")
      survey.sf <- sf::st_transform(survey.sf3, sf::st_crs(MAP$akland))[1:(nrow(survey.sf3) - 1), ]
    }
  }
  
  if (class(density.map)[1] == "sf") {
    density.sf <- sf::st_transform(density.map, sf::st_crs(MAP$akland))
    names(density.sf[1]) <- "density"
  }
  
  if (class(density.map)[1] == "RasterLayer") {
    density.dat0 <- as.data.frame(raster::xyFromCell(density.map, 1:raster::ncell(density.map)))
    vals <- raster::getValues(density.map)
    density.dat <- data.frame(lat = density.dat0$y, lon = density.dat0$x, density = vals)
    density.sf0 <- sf::st_as_sf(x = subset(density.dat, is.na(density) == FALSE), coords = c("lon", "lat"), crs = 3338)
    density.sf <- sf::st_transform(density.sf0, sf::st_crs(MAP$akland))
  }
  
  upper <- stats::quantile(density.sf$density, buffer, na.rm = TRUE)
  density.sf$density[density.sf$density > upper] <- upper
  
  densityplot <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = density.sf, ggplot2::aes(col = density), size = 0.05) + 
    ggplot2::geom_sf(data = survey.sf, fill = NA) + 
    ggplot2::geom_sf(data = MAP$akland, fill = "grey40") + 
    ggplot2::geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) + 
    ggplot2::geom_sf(data = MAP$bathymetry, color = "grey60", size = 0.25) + 
    ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust[1:2], ylim = MAP$plot.boundary$y + ext.adjust[3:4]) + 
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) + 
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) + 
    viridis::scale_color_viridis(option = col.palette, begin = col.palette.limits[1], end = col.palette.limits[2], na.value = NA, name = legend.title, labels = scales::comma_format(big.mark = ",")) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA), panel.background = ggplot2::element_rect(fill = NA, color = "black"), legend.key = ggplot2::element_rect(fill = NA, color = "grey30"), legend.position = legend.pos, legend.margin = ggplot2::margin(0, 0, 0, 0), axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12), legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12), plot.background = ggplot2::element_rect(fill = NA, color = NA))
  
  if (!is.na(title.name) && length(title.pos) == 2) {
    densityplot <- densityplot + ggplot2::geom_label(data = data.frame(x = title.pos[1], y = title.pos[2], label = title.name), ggplot2::aes(x = x, y = y, label = label, hjust = 0, vjust = 1), size = 5)
  }
  
  if (length(legend.pos) == 2) {
    densityplot <- densityplot + ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = barheight))
  }
  
  return(densityplot)
}
