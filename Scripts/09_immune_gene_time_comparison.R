# ==========================================
# Immune Gene Expression Comparison
# Time points: 3h, 48h and 96h
# Dataset: GSE43856
# ==========================================

library(ggplot2)

# Define immune-related genes
immune_genes <- c(
  "RUNX3", "CD247", "CD3D", "STAT4",
  "NFATC1", "PPP3CC", "IL2RB",
  "MAPK14", "HLA-DRA", "HLA-DMB"
)

# Read differential expression results
deg_3h <- read.csv("results/Annotated_DEG_3h.csv")
deg_48h <- read.csv("results/Annotated_DEG_48h.csv")
deg_96h <- read.csv("results/Annotated_DEG_96h.csv")

# Function to extract gene expression values
extract_immune <- function(data, time_point) {
  
  selected <- data[
    data$Symbol %in% immune_genes,
    c("Symbol", "logFC")
  ]
  
  selected <- selected[!is.na(selected$logFC), ]
  
  # Calculate median logFC for duplicate probes
  result <- aggregate(
    logFC ~ Symbol,
    data = selected,
    FUN = median
  )
  
  result$Time <- time_point
  
  return(result)
}

# Extract results for each time point
results_3h <- extract_immune(deg_3h, "3h")
results_48h <- extract_immune(deg_48h, "48h")
results_96h <- extract_immune(deg_96h, "96h")

# Combine results
immune_time_results <- rbind(
  results_3h,
  results_48h,
  results_96h
)

# Save comparison results
write.csv(
  immune_time_results,
  "results/Immune_gene_time_comparison.csv",
  row.names = FALSE
)

# Display results
print(immune_time_results)

# Create visualization
time_plot <- ggplot(
  immune_time_results,
  aes(
    x = Time,
    y = logFC,
    group = Symbol,
    color = Symbol
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Immune-related Gene Expression Across Time",
    subtitle = "GSE43856 | 3h, 48h and 96h After Exercise",
    x = "Time After Exercise",
    y = "Median Log Fold Change",
    color = "Gene"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Display plot
print(time_plot)

# Save visualization
ggsave(
  "figures/11_immune_gene_time_comparison.png",
  plot = time_plot,
  width = 12,
  height = 8,
  dpi = 300
)




immune_48h <- read.csv("results/Annotated_DEG_48h.csv")

immune_96h <- read.csv("results/Annotated_DEG_96h.csv")

immune_sig_48h <- immune_48h[
  immune_48h$Symbol %in% immune_genes &
    immune_48h$adj.P.Val < 0.05,
]

immune_sig_96h <- immune_96h[
  immune_96h$Symbol %in% immune_genes &
    immune_96h$adj.P.Val < 0.05,
]

cat("Significant immune-related probes at 48h:",
    nrow(immune_sig_48h), "\n")

cat("Significant immune-related probes at 96h:",
    nrow(immune_sig_96h), "\n")

print(immune_sig_48h[, c("Symbol", "logFC", "adj.P.Val")])

print(immune_sig_96h[, c("Symbol", "logFC", "adj.P.Val")])