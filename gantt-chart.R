################################################################################
# gantt-chart.R
# 
# Jeff Stevens
#
#
# Description
# -----------
#
# A Gantt-charting function for multi-plate assay workflows.
################################################################################


library(ggplot2)


gantt_chart <- function(timings, by = "plate", title = "Assay schedule Gantt chart") {
  # timings:  Data frame with columns named "Plate", "Resource", "Start", and
  # "Stop".  "Resource" should be an ordered factor.
  #
  # by: What to plot on the y-axis ("plate" for plate number, or "resource" for
  # the resource used).
  
  # Convert times to minutes
  timings <-
    timings %>%
    mutate( Start = Start / 60, Stop = Stop / 60 )
  
  if ( identical(by, "plate") ) {
    
    breaks <- order(unique(timings$Plate))
    labels <- as.character(as.integer(breaks))
    gridlines <- c(0.5, breaks + 0.5)
    
    p <- 
      ggplot(timings) +
        geom_rect( aes( xmin = Start, xmax = Stop,
                        ymin = Plate - 0.5, ymax = Plate + 0.5,
                        fill = Resource), color = "black" ) +
        scale_y_continuous( breaks = breaks, labels = labels, trans = "reverse",
                            minor_breaks = gridlines ) +
        xlab("\nTime (minutes)") + ylab("Plate\n") +
        ggtitle( paste0(title, "\nBroken out by plate", "\n") ) +
        theme_bw() +
        ## This suppresses the middle grid lines...
        theme( panel.grid.major.y = element_blank() )
    
  } else if ( identical(by, "resource") ) {
    
    breaks <- seq_along(levels(timings$Resource))
    labels <- levels(timings$Resource)
    gridlines <- c(0.5, breaks + 0.5)
    
    p <- 
      ggplot(timings) +
        geom_rect( aes( xmin = Start, xmax = Stop,
                        ymin = as.integer(Resource) - 0.5, ymax = as.integer(Resource) + 0.5,
                        fill = as.factor(Plate)), color = "black" ) +
        scale_y_continuous( breaks = breaks, labels = labels, trans = "reverse",
                            minor_breaks = gridlines) +
        scale_fill_discrete(name = "Plate") +
        xlab("\nTime (minutes)") + ylab("Resource\n") +
        ggtitle( paste0(title, "\nBroken out by resource", "\n") ) +
        theme_bw() +
        theme( panel.grid.major.y = element_blank() )
    
  } else {
    
    stop("'by' must either be 'plate' or 'resource'.")
  }
  
  
  return(p)
}