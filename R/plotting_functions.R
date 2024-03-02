#' @export
raster_occur_gp <-
  function(z1,
           z2 = NULL,
           z3 = NULL,
           coord,
           col_lim = list(c(-3, 3), c(0,1), c(-3,3)),
           titles = c("True Latent Proc.", "Obs. Binary Proc.", "Post. Mean. Latent  Proc.")) {
    if(!is.list(col_lim)){
      col_lim <- list(col_lim, col_lim, col_lim)
    }
    else if(is.null(col_lim)){
      col_lim <- list(range(z1), range(z2), range(z3))
    }
    else if(length(col_lim)!=3){
        stop("col_lim must either be a vector giving lower and upper bounds,
             or a list containing lower and upper bounds for each of z1,
             z2, z3.")
    }
    if(is.null(titles)){
      titles <- c("","","")
    }
    nloc <- nrow(coord)
    # First field
    p1 <-
      ggplot2::ggplot(data = data.frame(z = z1, lon = coord[, 1], lat = coord[, 2]),
             ggplot2::aes(x = lon, y = lat)) +
      ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = FALSE) +
      ggplot2::scale_fill_gradientn(colours = fields::tim.colors(n = nloc),
                                    limits = col_lim[[1]]) +
      ggplot2::ggtitle(titles)
    if(is.null(z2) & is.null(z3)) return(cowplot::plot_grid(p1, nrow = 1))

    # Second field
    p2 <-
      ggplot2::ggplot(data = data.frame(z = z2, lon = coord[, 1], lat = coord[, 2]),
             ggplot2::aes(x = lon, y = lat)) +
      ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = FALSE) +
      ggplot2::scale_fill_gradientn(colours = fields::tim.colors(n = nloc),
                           limits = col_lim[[2]]) +
      ggplot2::ggtitle(titles[2])
    if(is.null(z3)) return(cowplot::plot_grid(p1, p2, nrow = 1))

    # Third field
    p3 <-
      ggplot2::ggplot(data = data.frame(z = z3, lon = coord[, 1], lat = coord[, 2]),
             ggplot2::aes(x = lon, y = lat)) +
      ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = FALSE) +
      ggplot2::scale_fill_gradientn(colours = fields::tim.colors(n = nloc),
                                    limits = col_lim[[3]]) +
      ggplot2::ggtitle(titles[3])
    return(cowplot::plot_grid(p1,p2, p3, nrow = 1))
  }



#' @export
scatter_occur_gp <-
  function(z1,
           z2 = NULL,
           z3 = NULL,
           coord,
           col_lim = list(c(-3, 3), c(0,1), c(-3,3)),
           titles = c("True Latent Proc.", "Obs. Binary Proc.", "Post. Mean. Latent  Proc.")) {
    if(!is.list(col_lim)){
      col_lim <- list(col_lim, col_lim, col_lim)
    }
    else if(is.null(col_lim)){
      col_lim <- list(range(z1), range(z2), range(z3))
    }
    else if(length(col_lim)!=3){
        stop("col_lim must either be a vector giving lower and upper bounds,
             or a list containing lower and upper bounds for each of z1,
             z2, z3.")
    }
    nloc <- nrow(coord)
    # First Field
    p1 <-
      ggplot2::ggplot(data = data.frame(z = z1, lon = coord[, 1], lat = coord[, 2]),
                      ggplot2::aes(x = lon, y = lat, colour = z)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradientn(colours = fields::tim.colors(n = nloc), limits = col_lim[[1]]) +
      ggplot2::ggtitle(titles[1])
    if(is.null(z2) & is.null(z3)) return(cowplot::plot_grid(p1, nrow = 1))
    # Second Field
    p2 <-
      ggplot2::ggplot(data = data.frame(z = z2, lon = coord[, 1], lat = coord[, 2]),
                      ggplot2::aes(x = lon, y = lat,colour = z)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradientn(colours = fields::tim.colors(n = nloc), limits = col_lim[[2]]) +
      ggplot2::ggtitle(titles[2])
    if(is.null(z3)) return(cowplot::plot_grid(p1, p2, nrow = 1))

    # Third Field
    p3 <-
      ggplot2::ggplot(data = data.frame(z = z3, lon = coord[, 1], lat = coord[, 2]),
                      ggplot2::aes(x = lon, y = lat, colour = z)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradientn(colours = fields::tim.colors(n = nloc), limits = col_lim[[3]]) +
      ggplot2::ggtitle(titles[3])

    return(cowplot::plot_grid(p1,p2, p3, nrow = 1))
  }

#' @export
raster_world <- function(z, coord,
                         xlim = NULL,
                         ylim = NULL,
                         use_coord_lim = FALSE,
                         overlay_world_map = TRUE,
                         ncdf_coord_format = FALSE) {
  nloc <- nrow(coord)
  if (overlay_world_map) {
    world_map <- ggplot2::map_data("world")
  }
  if(ncdf_coord_format){
    coord[,1] <- (coord[, 1] + 180) %% 360 - 180
  }
  if(use_coord_lim){
    xlim <- range(coord[,1])
    ylim <- range(coord[,2])
  }
  p <- ggplot2::ggplot(data = data.frame(
    z = z,
    lon = coord[, 1],
    lat = coord[, 2]
  ),
  ggplot2::aes(x = lon, y = lat)) +
    ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = FALSE) +
    ggplot2::scale_fill_gradientn(colours = fields::tim.colors(n = nloc)) +
    ggplot2::geom_polygon(
      data = world_map,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = NA,
      color = "black",
      size = 0.5
    )  +
    ggplot2::theme(panel.background = ggplot2::element_blank())
    if(!is.null(xlim)){
      p <- p + ggplot2::xlim(xlim)
    }
    if(!is.null(ylim)){
      p <- p + ggplot2::ylim(ylim)
    }
  return(p)
}


#' Grid plot
#'
#' @param z vector or list of vectors to plot on map
#' @param coord nx2 matrix or list of nx2 matrices of coordinates for raster
#' or scatterplots
#' @param type one of "raster" or "point" or list of these specifying whether to
#' make a raster plot or scatterplot
#' @param col_lim upper and lower bounds for color limits of plot or list of limits
#' If null, will be inferred from data.
#' @param titles title or list of titles for each plot. If null, will be inferred
#' @param xlim xlimits for plot (if null, will be inferred from data)
#' @param ylim y-limits for plot (if null will be inferred from data)
#' @param polygon_df data.frame with lon and lat columns for overlaying polygons
#'  of,e.g. states
#'
#' @export
map_grid_plot <-  function(z,
                       coord,
                       type = "raster",
                       col_lim = NULL,
                       titles = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       polygon_df = NULL) {
  if(length(z) > 16){
    stop("Too many plots for a single grid.")
  }
  if(is.data.frame(coord)|(!is.list(coord))){
    coord <- list(coord)
  }
  if(!is.list(z)){
    z <- list(z)
  }
  if((length(coord) ==1)&(length(z) >1)){
    coord <- rep(coord, length(z))
  }
  if(!is.list(type)){
    type <- list(type)
  }
  if(length(type)==1){
    type <- rep(type, length(z))
  }
  if(is.null(col_lim)){
    col_lim <- lapply(z, range)
  }
  else if(!is.list(col_lim)){
    col_lim <- list(col_lim)
  }
  if(is.null(titles)){
    titles <- as.list(rep("", length(z)))
  }
  plot_list <- list()
  for(i in 1:length(z)){
    nloc <- nrow(coord[[i]])
    if(is.null(xlim)){
      xlim <- range(coord[[i]][,1])
    }
    if(is.null(ylim)){
      ylim <- range(coord[[i]][,2])
    }
    if (type[[i]] == "raster") {
      plot_list[[i]] <-
        ggplot2::ggplot(data = data.frame(z = z[[i]], lon = coord[[i]][, 1], lat = coord[[i]][, 2]),
                        ggplot2::aes(x = lon, y = lat)) +
        ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = FALSE) +
        ggplot2::scale_fill_gradientn(colours = fields::tim.colors(n = nloc),
                                      limits = col_lim[[i]]) +
        ggplot2::coord_quickmap(xlim = xlim,
                       ylim = ylim)
    }
    else{  # point type
      plot_list[[i]] <-
        ggplot2::ggplot(data = data.frame(z = z[[i]], lon = coord[[i]][, 1], lat = coord[[i]][, 2]),
                        ggplot2::aes(x = lon, y = lat)) +
        ggplot2::geom_point(ggplot2::aes(colour = z)) +
        ggplot2::scale_color_gradientn(colours = fields::tim.colors(n = nloc),
                                       limits = col_lim[[i]]) +
        ggplot2::coord_quickmap(xlim = xlim,
                       ylim = ylim)
    }
    plot_list[[i]] <- plot_list[[i]] +
      ggplot2::ggtitle(titles[[i]]) +
      ggplot2::theme(
        plot.title = element_text(size = 25,
                                  # family = "Arial Narrow",
                                  face = "plain"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        text = element_text(size = 10),
        legend.key.size = unit(1, "cm"),
        legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

    if(!is.null(polygon_df)){
      plot_list[[i]] <- plot_list[[i]] +
        ggplot2::geom_polygon(data = polygon_df,
                     ggplot2::aes(x = lon, y = lat, group = group),
                     fill = NA,
                     color = "black", size = 0.25)
    }
  }
  return(cowplot::plot_grid(plotlist = plot_list, nrow = ceiling(length(z)/3)))
}

