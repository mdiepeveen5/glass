recursiveCorPlot <- function(normalised_correlated_data, labels, font_scale , legend_scale , method="ward.D2", return_h_object = FALSE, caption=NULL) {
  col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                        "#4393C3", "#2166AC", "#053061"))
  
  
  # remove duplicate entries:
  plt <- normalised_correlated_data %>%
    tibble::rownames_to_column('__hugo_symbol__') %>%
    dplyr::filter(!duplicated(.data$`__hugo_symbol__`)) %>%
    tibble::column_to_rownames('__hugo_symbol__')
  
  
  #h <- hclust(dist(plt %>% as.matrix %>% t() %>%  scale %>% t() %>% as.matrix), method = method ) # Euclidean distance based clustering
  
  # determine correlation
  plt <- plt %>%
    as.matrix %>%
    t() %>%
    stats::cor()
  
  
  # find order by taking correlation of the correlation
  h <- stats::hclust(stats::as.dist(1 - stats::cor(plt)), method = method ) # recursive cor-based cluastering !!!
  #h <- stats::hclust( stats::as.dist(1 - plt) , method = method ) # regular cor-based clustering
  
  o <- h$labels[h$order] %>% rev()
  
  ph <- ggdendro::ggdendrogram(h, rotate = TRUE, theme_dendro = FALSE) +
    ggdendro::theme_dendro()
  
  # re-order to cor-cor clustering order and transform from data.frame into matrix
  plt <- plt %>%
    as.data.frame %>%
    dplyr::select(o) %>%
    t() %>%
    as.data.frame %>%
    dplyr::select(o) %>%
    t() %>%
    as.matrix
  
  
  # to test:
  # corrplot::corrplot(plt)
  
  
  
  # add x and y co-ordinates later on to melted table
  o.join <- data.frame(name = o, i = 1:length(o))
  
  plt.expanded2 <- reshape2::melt(plt) %>%
    dplyr::rename(y = .data$`Var1`) %>%
    dplyr::rename(x = .data$`Var2`) %>%
    dplyr::mutate(x = as.factor(.data$`x`)) %>%
    dplyr::mutate(y = as.factor(.data$`y`)) %>%
    dplyr::left_join(o.join %>% dplyr::rename(x.order = .data$`i`), by=c('x' = 'name'))%>%
    dplyr::left_join(o.join %>% dplyr::mutate(i = dplyr::n() - .data$i + 1  ) %>% dplyr::rename(y.order = .data$i), by=c('y' = 'name'))
  
  rm(o.join)
  
  
  
  p1 <- ggplot(
    plt.expanded2,
    aes(
      x = .data$x.order,
      y = .data$y.order,
      radius = ((abs(.data$value) * 0.7) + 0.3) / 2 - 0.05 ,
      # [0.3 , 0.8] + 0.2 smoothened from lwd/border
      fill = .data$value,
      col = .data$value,
      label = .data$x
    )
  ) +
    geom_tile( col="gray", fill="white", lwd=0.15) +
    scale_fill_gradientn( colours = col2(200), na.value = "grey50", limits = c(-1,1) , guide="none") + # guide = "colourbar",
    scale_color_gradientn( colours = col2(200), na.value = "grey50", limits = c(-1,1) , guide="none" ) +
    geom_circle(radius.fixed = T) + # from THIS repo
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme(legend.position = 'bottom',
          axis.text.y = element_text(size = font_scale, angle = 0, hjust = 1, vjust = 0.5), # used to be [3,6] reduce font size here, should become argument
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, color="gray80"),
          
          text = element_text(size=13),
          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank()
    ) +
    labs(y = NULL, x=NULL, main=NULL) +
    ggplot2::coord_fixed() +
    scale_y_continuous(name=NULL, breaks = length(o):1, labels = o)
  
  
  
  plt <- data.frame(gid = o, i = 1:length(o)) %>%
    dplyr::left_join(labels %>%
                       tibble::rownames_to_column('gid'), by=c('gid' = 'gid')) %>%
    reshape2::melt(id.vars = c('gid','i'))  %>%
    dplyr::mutate(variable = factor(.data$variable, levels = rev(colnames(labels))))
  
  
  p2 <- ggplot(plt , aes(x = .data$i ,
                         y = .data$variable ,
                         fill = .data$value,
                         label = .data$gid)) +
    geom_tile(col='white',lwd=0.15) +
    #scale_x_discrete(position = "bottom")  +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.text.y = element_text(size = font_scale, angle = 0, hjust = 1, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank()
    ) +
    guides(fill="none") +
    ggplot2::coord_fixed(ratio = legend_scale) + # used to be 2.75
    labs(x=NULL, y=NULL) +
    scale_fill_manual(values=c('TRUE'='red','FALSE'='gray98'))
  #scale_fill_manual(values=c('TRUE'='gray40','FALSE'='gray95'))
  
  
  #(p2 / p1) / ph
  
  #(p2 + plot_layout(guides = 'collect') ) /
  #(p1 + (ph))
  
  
  
  
  if(return_h_object) {
    return(h) # return clust object
  } else {
    layout <- '
A#
BC'
    
    return(
      patchwork::wrap_plots(
        A = p2,
        B = p1,
        C = (ph + plot_spacer () ),
        design = layout) +
        patchwork::plot_annotation(
          caption = caption
        )
    )
  }
}


scaleRadius <- function(data) {
  # radius -> [0,1], then to [0.01, 0.1]
  minr <- min(data$radius, na.rm = TRUE)
  maxr <- max(data$radius, na.rm = TRUE)
  if (minr == maxr) {
    if (minr > 0) {
      if (minr > 1) data$radius <- 1
      return (data)
    }
  }
  
  data$radius <- with(data, (radius - minr)/abs(maxr - minr))
  data$radius <- with(data, 0.01 + 0.1*radius)
  data
}

#' @rdname geom_circle
#' @export
GeomCircle <- ggplot2::ggproto("GeomCircle", ggplot2::Geom,
                               required_aes = c("x", "y", "radius", "radius.fixed"),
                               default_aes = ggplot2::aes(
                                 colour = "grey30", fill=NA, alpha=NA, linewidth=1, linetype="solid", radius=0.5 ),
                               
                               draw_key = function (data, params, size)
                               {
                                 grid::circleGrob(
                                   0.5, 0.5,
                                   r=0.35,
                                   gp = grid::gpar(
                                     col = scales::alpha(data$colour, data$alpha),
                                     fill = scales::alpha(data$fill, data$alpha),
                                     lty = data$linetype,
                                     lwd = data$linewidth,
                                     fontsize = data$radius*2
                                   )
                                 )
                               },
                               
                               
                               draw_panel = function(data, panel_scales, coord,  radius.fixed, na.rm = TRUE) {
                                 
                                 if(radius.fixed) {
                                   dx <- abs(panel_scales$x.range[2] - panel_scales$x.range[1])
                                   dy <- abs(panel_scales$y.range[2] - panel_scales$y.range[1])
                                   d <- min(dx, dy)
                                   
                                   coords <- coord$transform(data, panel_scales)
                                   
                                   coords$radius <- coords$radius / d # scale proportional to grid
                                   #coords <- scaleRadius(coords) # scaling to reset to [0, 1] unnecessary
                                 } else {
                                   coords <- coord$transform(data, panel_scales)
                                   coords <- scaleRadius(coords)
                                 }
                                 
                                 grid::circleGrob(
                                   x=coords$x, y=coords$y,
                                   r=coords$radius,
                                   gp = grid::gpar(
                                     col = alpha(coords$colour, coords$alpha),
                                     fill = alpha(coords$fill, coords$alpha),
                                     lty = coords$linetype,
                                     lwd = coords$linewidth,
                                     fontsize = coords$radius*2
                                   )
                                 )
                               }
)

#' Geom for drawing circles in the ggplot2 framework
#'
#' Circles are drawn with a specified radius centered at (x, y).
#' This geom is very much exploratory - we are using it for drawing edges for self references.
#' It is not explored for any more general use, so use with caution!
#' @inheritParams ggplot2::geom_point
#' @param radius numeric value giving the radius of the circle to be drawn (0-1 normalized scale)
#' @param radius.fixed Make the size of the radius fixed to grid coordinates instead of xlim and ylim (T/F)
#' @export
#' @examples
#' # circles are drawn centered at x and y
#' library(ggplot2)
#' data(mpg)
#' ggplot(mpg, aes(displ, hwy, radius=0.1)) +
#'   geom_circle(radius=0.1) + geom_point()
#' ggplot(mpg, aes(displ, hwy, radius=0.05)) +
#'   geom_circle(linetype=2, radius=0.05, alpha=0.5)
#' ggplot(mpg, aes(displ, hwy, radius=0.05)) +
#'   geom_circle(aes(linetype=factor(cyl)), radius=0.05, alpha=0.5)
#' df = data.frame(x = 0, y = 0)
#' ggplot(df, aes(x=x, y=y, radius=1)) + geom_point(cex=4) +
#'   geom_circle(radius=1, col="red", radius.fixed=TRUE) + xlim(-3,3) + ylim(-3,3)
#' ggplot(df, aes(x=x, y=y, radius=0.55)) + geom_point(cex=4) +
#'   geom_circle(radius=0.55, col="red", radius.fixed=FALSE) + xlim(-3,3) + ylim(-3,3)

geom_circle <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, radius.fixed = FALSE, ...) {
  
  ggplot2::layer(
    geom = GeomCircle, mapping = mapping,  data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, radius.fixed=radius.fixed, ...)
  )
}

primrecurRNA %>% dplyr::select("104059-002-001")
