# ==========================================
# KEGG Pathway Enrichment Analysis
# Dataset: GSE43856
# Time point: 3 hours after exercise
# ==========================================

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load annotated differential expression results
deg_3h <- read.csv("results/Annotated_DEG_3h.csv")

# Select genes with valid Entrez IDs
deg_valid <- deg_3h[
  !is.na(deg_3h$Entrez_Gene_ID) &
    deg_3h$Entrez_Gene_ID != "",
]

# Select significantly differentially expressed genes
sig_genes <- deg_valid[deg_valid$adj.P.Val < 0.05, ]

# Prepare gene lists
sig_entrez <- unique(as.character(sig_genes$Entrez_Gene_ID))

background_entrez <- unique(as.character(deg_valid$Entrez_Gene_ID))

# Perform KEGG enrichment analysis
kegg_result <- enrichKEGG(
  gene = sig_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20
)

# Convert Entrez IDs into gene symbols
kegg_result <- setReadable(
  kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# Convert results into a data frame
kegg_table <- as.data.frame(kegg_result)

# Display results
print(head(kegg_table, 10))

# Save results
write.csv(
  kegg_table,
  "results/KEGG_enrichment_3h.csv",
  row.names = FALSE
)

# Check number of enriched pathways
cat("Number of enriched KEGG pathways:", nrow(kegg_table))



# ==========================================
# KEGG Enrichment Visualization
# ==========================================

# Create output directory if it does not exist
dir.create("figures", showWarnings = FALSE)

# Generate dot plot
kegg_plot <- dotplot(
  kegg_result,
  showCategory = 8,
  font.size = 10
) +
  labs(
    title = "KEGG Pathway Enrichment",
    subtitle = "GSE43856 | 3 Hours After Exercise",
    x = "Gene Ratio",
    y = "KEGG Pathways"
  ) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    )
  )

# Display plot
print(kegg_plot)

# Save plot
ggsave(
  "figures/09_KEGG_enrichment_3h.png",
  plot = kegg_plot,
  width = 13,
  height = 9,
  dpi = 300
)



kegg_table$Description
kegg_table[, c("Description", "GeneRatio", "p.adjust", "geneID")]

immune_genes <- c(
  "RUNX3", "CD247", "CD3D", "STAT4",
  "NFATC1", "PPP3CC", "IL2RB",
  "MAPK14", "HLA-DRA", "HLA-DMB"
)

immune_results <- sig_genes[
  sig_genes$Symbol %in% immune_genes,
  c("Symbol", "logFC", "P.Value", "adj.P.Val")
]

immune_results <- immune_results[
  order(immune_results$adj.P.Val),
]

print(immune_results)


immune_unique <- immune_results[
  order(immune_results$adj.P.Val),
]

immune_unique <- immune_unique[
  !duplicated(immune_unique$Symbol),
]

print(immune_unique)


write.csv(
  immune_unique,
  "results/Immune_related_genes_3h.csv",
  row.names = FALSE
)


# ==========================================
# Immune-related Gene Expression Visualization
# ==========================================

library(ggplot2)

# Arrange genes according to logFC
immune_unique$Symbol <- factor(
  immune_unique$Symbol,
  levels = immune_unique$Symbol[
    order(immune_unique$logFC)
  ]
)

# Create bar plot
immune_plot <- ggplot(
  immune_unique,
  aes(x = Symbol, y = logFC, fill = logFC > 0)
) +
  geom_col() +
  coord_flip() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
    labels = c("Downregulated", "Upregulated")
  ) +
  labs(
    title = "Immune-related Gene Expression Changes",
    subtitle = "GSE43856 | 3 Hours After Exercise",
    x = "Gene",
    y = "Log Fold Change",
    fill = "Expression"
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
print(immune_plot)

# Save plot
ggsave(
  "figures/10_immune_gene_expression_3h.png",
  plot = immune_plot,
  width = 10,
  height = 7,
  dpi = 300
)