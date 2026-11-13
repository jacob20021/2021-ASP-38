if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("GEOquery")
BiocManager::install("limma")

install.packages(c(
  "tidyverse",
  "ggplot2",
  "pheatmap"
))
library(GEOquery)
library(limma)
library(tidyverse)
library(ggplot2)
library(pheatmap)

install.packages("farver", dependencies = TRUE)
library(farver)
library(ggplot2)

install.packages(
  "farver",
  repos = "https://cloud.r-project.org",
  dependencies = FALSE,
  type = "binary"
)

library(farver)

library(ggplot2)
library(tidyverse)
library(pheatmap)

gse <- getGEO("GSE43856", GSEMatrix = TRUE)

raw_file <- "C:/Users/anopj/Downloads/GSE43856_RAW.tar"

file.exists(raw_file)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

untar(raw_file, exdir = "data/raw")

list.files("data/raw", recursive = TRUE)
list.files("data/raw", recursive = TRUE, full.names = TRUE)

getwd()
list.files("data/raw", recursive = TRUE)
list.files("data/raw", recursive = TRUE, full.names = TRUE)

con <- gzfile("data/raw/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz", "rt")

readLines(con, n = 10)

close(con)

con <- gzfile("data/raw/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz", "rt")

readLines(con, n = 10)

close(con)

file.info(list.files("data/raw", full.names = TRUE))$size
file1 <- "data/raw/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz"

file2 <- "data/raw/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"

con1 <- gzfile(file1, "rt")
lines1 <- readLines(con1, n = 10)
close(con1)

con2 <- gzfile(file2, "rt")
lines2 <- readLines(con2, n = 10)
close(con2)

cat("FIRST FILE:\n")
cat(lines1, sep = "\n")

cat("\n\nSECOND FILE:\n")
cat(lines2, sep = "\n")

file1 <- "data/raw/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz"

con <- gzfile(file1, "rt")

lines1 <- readLines(con, n = 10)

close(con)

cat(lines1, sep = "\n")

list.files("data/raw", recursive = TRUE)

list.files("data/raw", recursive = TRUE)


file.exists("data/raw/GSE43856_series_matrix.txt.gz")

library(GEOquery)

gset <- getGEO(
  filename = "data/raw/GSE43856_series_matrix.txt.gz"
)

class(gset)

dim(exprs(gset))
dim(pData(gset))
head(pData(gset))

> 
> gset <- getGEO(
+   filename = "data/raw/GSE43856_series_matrix.txt.gz"
+ )
> library(GEOquery)
> 
> gset <- getGEO(
+   filename = "data/raw/GSE43856_series_matrix.txt.gz"
+ )
 # Extract sample information
 metadata <- pData(gset)
 
 # Count samples by tissue
 table(metadata$tissue.ch1)

 
 
 # Count samples by time point
 table(metadata$`time point:ch1`)
 
 
 # Count samples by tissue
 table(metadata$`tissue:ch1`)
 
 # Count samples by tissue and time point
 table(metadata$`tissue:ch1`, metadata$`time point:ch1`)
 
 
 # Select only blood PMN samples
 pmn_samples <- metadata$`tissue:ch1` == "blood PMN"
 
 # Extract neutrophil expression data
 expr_pmn <- exprs(gset)[, pmn_samples]
 
 # Extract neutrophil sample information
 meta_pmn <- metadata[pmn_samples, ]
 
 # Check dimensions
 dim(expr_pmn)
 dim(meta_pmn)
 
 
 
 # Select blood PMN samples
 pmn_samples <- metadata$`tissue:ch1` == "blood PMN"
 
 # Extract gene expression data
 expr_pmn <- exprs(gset)[, pmn_samples]
 
 # Extract sample information
 meta_pmn <- metadata[pmn_samples, ]
 
 # Check dimensions
 dim(expr_pmn)
 dim(meta_pmn)
 
 
 
 # Check first 6 gene probes and first 5 samples
 expr_pmn[1:6, 1:5]
 
 
 # Summary of expression values
 summary(as.vector(expr_pmn))
 
 sum(is.na(expr_pmn))
 
 
 expr_pmn[1:6, 1:5]
 
 summary(as.vector(expr_pmn))
 
 
 
 sum(is.na(expr_pmn))
 
 
 
 class(expr_pmn)
 is.numeric(expr_pmn)
 
 
 
 # Create a boxplot for all 32 samples
 boxplot(
   expr_pmn,
   main = "Gene Expression Distribution Across Neutrophil Samples",
   xlab = "Samples",
   ylab = "Log2 Expression",
   las = 2,
   outline = FALSE,
   cex.axis = 0.6
 )
 dir.create("figures", showWarnings = FALSE)
 
 dir.create(
   "C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures",
   showWarnings = FALSE
 )
 
 png(
   "C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures/01_gene_expression_distribution.png",
   width = 1600,
   height = 900,
   res = 150
 )
 
 boxplot(
   expr_pmn,
   main = "Gene Expression Distribution Across Neutrophil Samples",
   xlab = "Samples",
   ylab = "Expression Level",
   las = 2,
   outline = FALSE,
   cex.axis = 0.6
 )
 
 dev.off()
 file.exists("C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures/01_gene_expression_distribution.png")
 png(
   "C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures/01_gene_expression_distribution.png",
   width = 1800,
   height = 1200,
   res = 150
 )
 
 par(mar = c(10, 5, 4, 2))
 
 boxplot(
   expr_pmn,
   main = "Gene Expression Distribution Across Neutrophil Samples",
   xlab = "",
   ylab = "Log2 Expression",
   las = 2,
   outline = FALSE,
   cex.axis = 0.65
 )
 
 dev.off()
 
 figure_path <- "C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures"
 
 dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
 
 dir.exists(figure_path)
 
 
 
 png(
   filename = file.path(figure_path, "01_gene_expression_distribution_v2.png"),
   width = 1800,
   height = 1200,
   res = 150
 )
 
 par(mar = c(10, 5, 4, 2))
 
 boxplot(
   expr_pmn,
   main = "Gene Expression Distribution Across Neutrophil Samples",
   xlab = "",
   ylab = "Log2 Expression",
   las = 2,
   outline = FALSE,
   cex.axis = 0.65
 )
 
 dev.off()
 
 file.exists(
   "C:/Users/anopj/OneDrive/Documents/Bioinformatics/figures/01_gene_expression_distribution_v2.png"
 )