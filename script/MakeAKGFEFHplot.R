MakeAKGFEFHplot = function (region, efh.map, drop = 10^8, survey.area = "default", 
          ext.adjust = "default", legend.pos = "default", legend.labels = c("95%", 
                                                                            "75%", "50%", "25%"), legend.title = NA, col.palette = "viridis", 
          col.palette.limits = c(0, 1), title.name = NA, title.pos = "default") 
{
  region <- tolower(region)
  if (region %in% c("ebs", "bs.all", "sebs", "bs.south", "ecs", 
                    "ebs.ecs", "ai", "ai.west", "ai.central", "ai.east", 
                    "goa", "goa.west", "goa.east")) {
    MAP <- akgfmaps::get_base_layers(select.region = region, 
                                     set.crs = "auto")
  }
  else {
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
        legend.pos <- c(0.07, 0.28)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", 
                        "ecs", "ebs.ecs")) {
        legend.pos <- c(0.12, 0.18)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        legend.pos <- c(0.08, 0.51)
      }
    }
    else {
      legend.pos <- NA
    }
  }
  if (length(title.pos) < 2) {
    if (title.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", 
                        "ai.east")) {
        title.pos <- c(-9e+05, 490000)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", 
                        "ecs", "ebs.ecs")) {
        title.pos <- c(-550000, 8e+05)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        title.pos <- c(-1300000, 630000)
      }
    }
    else {
      title.pos <- NA
    }
  }
  if (is.character(survey.area) && survey.area[1] == "default") {
    survey.sf <- MAP$survey.area
  }
  else {
    if (class(survey.area)[1] == "sf") {
      survey.sf <- sf::st_transform(survey.area, sf::st_crs(MAP$akland))[1:(nrow(survey.area) - 
                                                                              1), ]
    }
    if (class(survey.area)[1] == "RasterLayer") {
      survey.sf1 <- stars::st_as_stars(is.na(survey.area))
      survey.sf2 <- sf::st_as_sf(survey.sf1, merge = TRUE)
      survey.sf3 <- sf::st_cast(survey.sf2, "POLYGON")
      survey.sf <- sf::st_transform(survey.sf3, sf::st_crs(MAP$akland))[1:(nrow(survey.sf3) - 
                                                                             1), ]
    }
  }
  efh.vals <- raster::getValues(efh.map)
  efh.vals[efh.vals == 1] <- NA
  efhpoly0 <- stars::st_as_stars(efh.map)
  efhpoly <- sf::st_as_sf(efhpoly0, merge = TRUE)
  efhpoly2 <- efhpoly[efhpoly$layer != 1, ]
  efh.dummy.raster <- raster::raster(efh.map)
  efh.vals2 <- is.na(efh.vals) == F
  efh.dummy.raster <- raster::setValues(efh.dummy.raster, values = efh.vals2)
  efhdummy0 <- stars::st_as_stars(efh.dummy.raster)
  efhdummy <- sf::st_cast(sf::st_as_sf(efhdummy0, merge = TRUE))
  efhdummy2 <- sf::st_transform(efhdummy, sf::st_crs(MAP$akland))
  efhdummy.poly <- sf::st_cast(efhdummy2, "POLYGON")
  areas <- sf::st_area(efhdummy.poly)
  outside <- order(areas, decreasing = T)[1]
  toosmall <- which(as.numeric(areas) < drop)
  efh.x <- efhdummy2$layer[-c(outside, toosmall)]
  efh.y <- efhdummy2$geometry[-c(outside, toosmall)]
  efhdummy3 <- sf::st_sf(efh.x, efh.y)
  efhplot <- ggplot2::ggplot() + ggplot2::geom_sf(data = survey.sf, 
                                                  fill = "grey95") + ggplot2::geom_sf(data = efhpoly2, 
                                                                                      ggplot2::aes(fill = as.factor(layer)), col = NA) + ggplot2::geom_sf(data = efhdummy3, 
                                                                                                                                                          fill = NA, size = 0.3) + ggplot2::geom_sf(data = MAP$akland, 
                                                                                                                                                                                                    fill = "grey40") + ggplot2::geom_sf(data = MAP$graticule, 
                                                                                                                                                                                                                                        color = "grey70", alpha = 0.5) + ggplot2::geom_sf(data = MAP$bathymetry, 
                                                                                                                                                                                                                                                                                          color = "grey60", size = 0.25)
  efhplot <- efhplot + ggplot2::coord_sf(xlim = MAP$plot.boundary$x + 
                                           ext.adjust[1:2], ylim = MAP$plot.boundary$y + ext.adjust[3:4]) + 
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) + 
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) + 
    viridis::scale_fill_viridis(discrete = T, name = legend.title, 
                                labels = legend.labels) + ggplot2::theme_bw() + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", 
                                                                                                                                    fill = NA), panel.background = ggplot2::element_rect(fill = NA, 
                                                                                                                                                                                         color = "black"), legend.key = ggplot2::element_rect(fill = NA, 
                                                                                                                                                                                                                                              color = "grey30"), legend.position = legend.pos, axis.title = ggplot2::element_blank(), 
                                                                                               axis.text = ggplot2::element_text(size = 12), legend.text = ggplot2::element_text(size = 12), 
                                                                                               legend.title = ggplot2::element_text(size = 12), plot.background = ggplot2::element_rect(fill = NA, 
                                                                                                                                                                                        color = NA))
  if (is.na(title.name) == F && length(title.pos) == 2) {
    efhplot <- efhplot + ggplot2::geom_label(data = data.frame(x = title.pos[1], 
                                                               y = title.pos[2], label = title.name), ggplot2::aes(x = x, 
                                                                                                                   y = y, label = label, hjust = 0, vjust = 1), size = 5)
  }
  if (length(legend.pos) == 2) {
    efhplot <- efhplot + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, 
                                                                                           shape = 15)))
  }
  return(efhplot)
}