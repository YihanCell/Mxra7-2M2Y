# Install and load necessary packages
library(ggplot2)
library(gridExtra)

# Prepare the data for WT2m, X2m, WT2y, and X2y
prop.2m <- data.frame(
  Sample = rep(c("WT2m", "X2m"), each = 5),
  Gene = c("B.FrF", "preB.FrD", "B.FrE", "preB.FrC", "proB.FrA"),
  Proportion = c(0.11392620, 0.06104404, 0.15592586, 0.23329366, 0.43581024,
                 0.17557932, 0.21372549, 0.18591800, 0.12976827, 0.29500891)
)

prop.2y <- data.frame(
  Sample = rep(c("WT2y", "X2y"), each = 5),
  Gene = c("preB.FrD", "B.FrF", "B.FrE", "preB.FrC", "proB.FrA"),
  Proportion = c(0.19379586, 0.54769847, 0.10273516, 0.11440961, 0.04136091,
                 0.16183507, 0.44610189, 0.11117366, 0.16549395, 0.11539544)
)

# Create a function to draw ring charts
draw_ring_chart <- function(data, title, colors) {
  ggplot(data, aes(x = 1, y = Proportion, fill = Gene)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = colors) +  # Set the custom color palette
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Specify the custom color palette (viridis)
viridis_colors <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")

# Draw the four ring charts
chart_wt2m <- draw_ring_chart(prop.2m[prop.2m$Sample == "WT2m", ], "WT2m", viridis_colors)
chart_x2m <- draw_ring_chart(prop.2m[prop.2m$Sample == "X2m", ], "X2m", viridis_colors)
chart_wt2y <- draw_ring_chart(prop.2y[prop.2y$Sample == "WT2y", ], "WT2y", viridis_colors)
chart_x2y <- draw_ring_chart(prop.2y[prop.2y$Sample == "X2y", ], "X2y", viridis_colors)

# Arrange all four charts in one picture
dev.new()
grid.arrange(chart_wt2m, chart_x2m, chart_wt2y, chart_x2y, ncol = 2)
dev.off()
