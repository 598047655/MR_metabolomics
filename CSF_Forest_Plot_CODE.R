# Load necessary library
library(forestploter)

# Clear all existing variables
rm(list = ls())

# Import the CSV file (ensure 'me.csv' is in the working directory)
df <- read.csv("me.csv")

# Add space column for formatting
df$' ' <- paste(rep(" ", 45), collapse = " ")

# Define the forest plot theme
tm <- forest_theme(
  base_size = 14,
  ci_pch = 16,
  ci_col = "#4575b4",
  ci_lty = 1,
  ci_lwd = 2,
  ci_Theight = 0.2,
  refline_lwd = 1.5,
  refline_lty = "dashed",
  refline_col = "grey20",
  footnote_cex = 0.8,
  footnote_fontface = "italic",
  footnote_col = "grey30",
  core = list(bg_params = list(fill = c("#FFFFFF", "#E6E6E6"), col = NA),
              fg_params = list(fontface = c(rep(c(2, 3, 3, 3, 3), 9)))),
  colhead = list(fg_params = list(col = "navyblue", fontface = 2))
)

# Create the forest plot
p <- forest(
  data = df[, c(1:3, 8, 4)],
  lower = df$low,
  upper = df$high,
  est = df$mean,
  ci_column = 4,
  sizes = df$mean,
  ref_line = 1,
  xlim = c(0.7, 1.8),
  ticks_at = c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8),
  theme = tm
)

# Print the plot to the console
print(p)

# Save the plot as a PNG file
png("forestploter.png", height = 50, width = 50, res = 300, units = "cm")
print(p)  # Ensure the plot is printed to the device
dev.off()
