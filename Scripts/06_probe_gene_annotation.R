
# 06_PROBE_GENE_ANNOTATION.R
# Probe-to-Gene Annotation for GSE43856

# 1. Set project directories

project_path <- "C:/Users/anopj/OneDrive/Documents/Bioinformatics"

data_path <- file.path(project_path, "data/raw")

results_path <- file.path(project_path, "results")

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)


# 2. Read Illumina annotation files

annotation_R1 <- read.delim(
  gzfile(file.path(
    data_path,
    "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz"
  )),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)

annotation_R2 <- read.delim(
  gzfile(file.path(
    data_path,
    "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"
  )),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)


# 3. Check annotation information

dim(annotation_R1)

dim(annotation_R2)

colnames(annotation_R1)

colnames(annotation_R2)


# 4. Check the first few rows

head(annotation_R1)

head(annotation_R2)


# 5. Identify common probe identifiers

intersect(
  colnames(annotation_R1),
  colnames(annotation_R2)
)


# 6. Check whether annotation files contain probe IDs

head(annotation_R1[[1]])

head(annotation_R2[[1]])


# Inspect the first 40 lines of annotation R1

con <- gzfile(file.path(
  data_path,
  "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz"
))

lines_R1 <- readLines(con, n = 40)
close(con)

cat(lines_R1, sep = "\n")


# Inspect the first 40 lines of annotation R2

con <- gzfile(file.path(
  data_path,
  "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"
))

lines_R2 <- readLines(con, n = 40)
close(con)

cat(lines_R2, sep = "\n")

#
annotation_R2 <- read.delim(
  gzfile(file.path(data_path,
                   "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz")),
  skip = 8,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(annotation_R2)
colnames(annotation_R2)
head(annotation_R2[, c("Search_Key", "Symbol", "Entrez_Gene_ID")])


# Check the first 10 non-control probes
human_probes <- annotation_R2[
  annotation_R2$Species == "Homo sapiens",
]

dim(human_probes)

head(
  human_probes[, c("Search_Key", "Symbol", "Entrez_Gene_ID")]
)

# Check the first 10 probe IDs in the DE results
head(rownames(results_3h), 10)

# Check the first 10 annotation probe IDs
head(human_probes$Search_Key, 10)

# Check annotation probe IDs
head(human_probes$Probe_Id, 10)

# Check whether DE probe IDs match annotation Probe_Id
sum(rownames(results_3h) %in% human_probes$Probe_Id)

# Check total DE probes
nrow(results_3h)




# Select required annotation columns
annotation_selected <- human_probes[, c(
  "Probe_Id",
  "Symbol",
  "Entrez_Gene_ID",
  "Definition"
)]

# Function to add gene annotations
annotate_results <- function(results) {
  
  results$Probe_Id <- rownames(results)
  
  annotated <- merge(
    results,
    annotation_selected,
    by = "Probe_Id",
    all.x = TRUE,
    sort = FALSE
  )
  
  return(annotated)
}

# Annotate all comparisons
annotated_3h <- annotate_results(results_3h)
annotated_48h <- annotate_results(results_48h)
annotated_96h <- annotate_results(results_96h)

# Save annotated results
write.csv(annotated_3h,
          file.path(results_path, "Annotated_DEG_3h.csv"),
          row.names = FALSE)

write.csv(annotated_48h,
          file.path(results_path, "Annotated_DEG_48h.csv"),
          row.names = FALSE)

write.csv(annotated_96h,
          file.path(results_path, "Annotated_DEG_96h.csv"),
          row.names = FALSE)



head(annotated_3h[, c(
  "Probe_Id",
  "Symbol",
  "Entrez_Gene_ID",
  "logFC",
  "P.Value",
  "adj.P.Val"
)])

sum(!is.na(annotated_3h$Symbol))


# Candidate immune and inflammatory genes
immune_genes <- c("IL6", "IL1B", "TNF", "CXCL8", "TLR4")

# Find these genes in the annotated results
immune_results <- annotated_3h[
  annotated_3h$Symbol %in% immune_genes,
]

immune_results[, c(
  "Probe_Id",
  "Symbol",
  "logFC",
  "P.Value",
  "adj.P.Val"
)]




# Check immune and inflammatory genes at 48h and 96h

immune_genes <- c("IL6", "IL1B", "TNF", "CXCL8", "TLR4")

immune_results_48h <- annotated_48h[
  annotated_48h$Symbol %in% immune_genes,
]

immune_results_96h <- annotated_96h[
  annotated_96h$Symbol %in% immune_genes,
]

# Display 48-hour results
immune_results_48h[, c(
  "Probe_Id", "Symbol", "logFC", "P.Value", "adj.P.Val"
)]

# Display 96-hour results
immune_results_96h[, c(
  "Probe_Id", "Symbol", "logFC", "P.Value", "adj.P.Val"
)]