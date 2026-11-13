# Match metadata with expression samples
meta_pmn <- meta_pmn[colnames(expr_pmn), ]

# Select the 500 most variable genes
gene_variance <- apply(expr_pmn, 1, var)

top_genes <- order(gene_variance, decreasing = TRUE)[1:500]

# Perform PCA
pca_result <- prcomp(
  t(expr_pmn[top_genes, ]),
  center = TRUE,
  scale. = FALSE
)

# Calculate percentage of variance
pca_variance <- round(
  100 * pca_result$sdev^2 / sum(pca_result$sdev^2),
  2
)

# Save PCA plot
png(
  "C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures/02_PCA_neutrophil_samples.png",
  width = 1600,
  height = 1200,
  res = 150
)

plot(
  pca_result$x[, 1],
  pca_result$x[, 2],
  col = as.numeric(as.factor(meta_pmn$`time point:ch1`)),
  pch = 19,
  xlab = paste("PC1 (", pca_variance[1], "%)", sep = ""),
  ylab = paste("PC2 (", pca_variance[2], "%)", sep = ""),
  main = "PCA of Neutrophil Gene Expression"
)

legend(
  "topright",
  legend = levels(as.factor(meta_pmn$`time point:ch1`)),
  col = 1:length(levels(as.factor(meta_pmn$`time point:ch1`))),
  pch = 19
)

dev.off()

