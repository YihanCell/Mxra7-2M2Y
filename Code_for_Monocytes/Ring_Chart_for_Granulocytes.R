# Install and load necessary packages
library(ggplot2)
library(gridExtra)

# Prepare the data for WT2m, X2m, WT2y, and X2y
prop.2y <- data.frame(
  Sample = rep(c("WT2y", "X2y"), each = 2),
  Gene = c("Neutrophils", "Basophils"),
  Proportion = c(0.97981940, 0.02018060,0.96599667,0.03400333)
)

prop.2m <- data.frame(
  Sample = rep(c("WT2m", "X2m"), each = 2),
  Gene = c("Neutrophils", "Basophils"),
  Proportion = c(0.96292095, 0.03707905, 0.97584396, 0.02415604)
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
viridis_colors <- c("#1f78b4", "#33a02c")

# Draw the four ring charts
chart_wt2m <- draw_ring_chart(prop.2m[prop.2m$Sample == "WT2m", ], "WT2m", viridis_colors)
chart_x2m <- draw_ring_chart(prop.2m[prop.2m$Sample == "X2m", ], "X2m", viridis_colors)
chart_wt2y <- draw_ring_chart(prop.2y[prop.2y$Sample == "WT2y", ], "WT2y", viridis_colors)
chart_x2y <- draw_ring_chart(prop.2y[prop.2y$Sample == "X2y", ], "X2y", viridis_colors)

# Arrange all four charts in one picture
dev.new()
grid.arrange(chart_wt2m, chart_x2m, chart_wt2y, chart_x2y, ncol = 2)
dev.off()


