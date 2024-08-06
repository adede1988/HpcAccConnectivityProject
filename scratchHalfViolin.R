# Defining the geom_flat_violin function ----
# Note: the below code modifies the
# existing github page by removing a parenthesis in line 50
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}
#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 1.5
              )
          },
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          draw_key = draw_key_polygon,
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "black", size = 2,
            alpha = 1, linetype = "solid"
          ),
          required_aes = c("x", "y")
  )


# Custom position function to combine jitter and nudge
position_jitternudge <- function(jitter.width = 0.2, jitter.height = 0, nudge.x = 0, nudge.y = 0) {
  ggproto(NULL, PositionJitternudge, 
          jitter.width = jitter.width, 
          jitter.height = jitter.height, 
          nudge.x = nudge.x, 
          nudge.y = nudge.y)
}

PositionJitternudge <- ggproto("PositionJitternudge", Position,
                               jitter.width = NULL,
                               jitter.height = NULL,
                               nudge.x = NULL,
                               nudge.y = NULL,
                               
                               setup_params = function(self, data) {
                                 list(
                                   jitter.width = self$jitter.width,
                                   jitter.height = self$jitter.height,
                                   nudge.x = self$nudge.x,
                                   nudge.y = self$nudge.y
                                 )
                               },
                               
                               compute_layer = function(data, params, panel) {
                                 data <- transform_position(data, params$jitter.width, params$jitter.height, "jitter")
                                 transform_position(data, params$nudge.x, params$nudge.y, "nudge")
                               }
)