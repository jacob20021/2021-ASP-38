deg_3h <- read.csv("results/Annotated_DEG_3h.csv")

names(deg_3h)
head(deg_3h)

# Total number of annotated genes
nrow(deg_3h)

# Number of significant genes
sum(deg_3h$adj.P.Val < 0.05, na.rm = TRUE)

# Significantly upregulated genes
sum(deg_3h$adj.P.Val < 0.05 & deg_3h$logFC > 0, na.rm = TRUE)

# Significantly downregulated genes
sum(deg_3h$adj.P.Val < 0.05 & deg_3h$logFC < 0, na.rm = TRUE)

# Check missing Entrez IDs
sum(is.na(deg_3h$Entrez_Gene_ID))

#--------------------
# Remove rows with missing Entrez IDs
deg_valid <- deg_3h[
  !is.na(deg_3h$Entrez_Gene_ID) &
    deg_3h$Entrez_Gene_ID != "",
]

# Significant genes
sig_genes <- deg_valid[deg_valid$adj.P.Val < 0.05, ]

# Unique significant Entrez IDs
sig_entrez <- unique(as.character(sig_genes$Entrez_Gene_ID))

# Background gene set (all tested genes)
background_entrez <- unique(as.character(deg_valid$Entrez_Gene_ID))

# Upregulated genes
up_entrez <- unique(as.character(
  sig_genes$Entrez_Gene_ID[sig_genes$logFC > 0]
))

# Downregulated genes
down_entrez <- unique(as.character(
  sig_genes$Entrez_Gene_ID[sig_genes$logFC < 0]
))

# Display counts
cat("Unique significant genes:", length(sig_entrez), "\n")
cat("Background genes:", length(background_entrez), "\n")
cat("Upregulated genes:", length(up_entrez), "\n")
cat("Downregulated genes:", length(down_entrez), "\n")




if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))




library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)




# ==========================================
# GO Biological Process Enrichment Analysis
# ==========================================

go_bp <- enrichGO(
  gene = sig_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# Convert results into a data frame
go_results <- as.data.frame(go_bp)

# Display the first 10 results
head(go_results, 10)

# Save GO enrichment results
write.csv(
  go_results,
  "results/GO_BP_enrichment_3h.csv",
  row.names = FALSE
)

# Plot GO enrichment results
if (nrow(go_results) > 0) {
  
  go_plot <- dotplot(
    go_bp,
    showCategory = 15,
    title = "GO Biological Process Enrichment - 3h"
  )
  
  print(go_plot)
  
  ggsave(
    "figures/07_GO_BP_enrichment_3h.png",
    plot = go_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
  
} else {
  message("No significantly enriched GO terms were found.")
}

install.packages("labeling")

# ==========================================
# GO Biological Process Enrichment Plot
# ==========================================

library(ggplot2)
library(enrichplot)

if (nrow(go_results) > 0) {
  
  # Create GO enrichment dot plot
  go_plot <- dotplot(
    go_bp,
    showCategory = 10,
    font.size = 10,
    label_format = 40
  ) +
    labs(
      title = "GO Biological Process Enrichment",
      subtitle = "GSE43856 | 3 Hours After Exercise",
      x = "Gene Ratio",
      y = "GO Biological Processes",
      color = "Adjusted P-value",
      size = "Gene Count"
    ) +
    theme(
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      
      plot.subtitle = element_text(
        size = 9,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      plot.margin = margin(8, 12, 8, 12)
    )
  
  # Display plot
  print(go_plot)
  
  # Save high-resolution plot
  ggsave(
    filename = "figures/07_GO_BP_enrichment_3h.png",
    plot = go_plot,
    width = 13,
    height = 9,
    units = "in",
    dpi = 300
  )
  
  cat("GO enrichment plot saved successfully!\n")
  
} else {
  
  cat("No significantly enriched GO terms were found.\n")
  
}


nrow(go_results)
head(go_results, 10)



write.csv(
  go_results,
  "results/GO_BP_enrichment_3h.csv",
  row.names = FALSE
)

# GO enrichment for upregulated genes
go_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# GO enrichment for downregulated genes
go_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# Save results
write.csv(
  as.data.frame(go_up),
  "results/GO_BP_upregulated_3h.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(go_down),
  "results/GO_BP_downregulated_3h.csv",
  row.names = FALSE
)

cat("Upregulated GO terms:", nrow(as.data.frame(go_up)), "\n")
cat("Downregulated GO terms:", nrow(as.data.frame(go_down)), "\n")

head(as.data.frame(go_down), 10)


# Plot downregulated GO enrichment

go_down_plot <- dotplot(
  go_down,
  showCategory = 10,
  font.size = 10,
  label_format = 40
) +
  labs(
    title = "GO Enrichment of Downregulated Genes",
    subtitle = "GSE43856 | 3 Hours After Exercise",
    x = "Gene Ratio",
    y = "GO Biological Processes"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

print(go_down_plot)

ggsave(
  "figures/08_GO_BP_downregulated_3h.png",
  plot = go_down_plot,
  width = 13,
  height = 9,
  dpi = 300
)