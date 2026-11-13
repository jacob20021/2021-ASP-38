# Create a correlation colour palette
cor_colors <- colorRampPalette(
  c("blue", "white", "red")
)(100)

# Save corrected heatmap
png(
  "C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures/03_sample_correlation_heatmap.png",
  width = 1800,
  height = 1600,
  res = 150
)

heatmap(
  sample_cor,
  scale = "none",
  col = cor_colors,
  zlim = c(-1, 1),
  symm = TRUE,
  main = "Sample Correlation Heatmap"
)

dev.off()