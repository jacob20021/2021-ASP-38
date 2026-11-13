# Load limma
library(limma)

# Prepare time point information
time_point <- sub(
  "^[0-3]; ",
  "",
  meta_pmn$`time point:ch1`
)

time_point <- factor(
  time_point,
  levels = c(
    "baseline pre-EXTRI",
    "3 hr post-EXTRI",
    "48 hr post-EXTRI",
    "96 hr post-EXTRI"
  )
)

# Participant information
subject_id <- factor(meta_pmn$`subject id:ch1`)

# Create design matrix
design <- model.matrix(~ 0 + time_point)

colnames(design) <- c(
  "Baseline",
  "H3",
  "H48",
  "H96"
)

# Account for repeated measurements
corfit <- duplicateCorrelation(
  expr_pmn,
  design,
  block = subject_id
)

# Fit linear model
fit <- lmFit(
  expr_pmn,
  design,
  block = subject_id,
  correlation = corfit$consensus
)

# Define comparisons against baseline
contrast_matrix <- makeContrasts(
  H3_vs_Baseline = H3 - Baseline,
  H48_vs_Baseline = H48 - Baseline,
  H96_vs_Baseline = H96 - Baseline,
  levels = design
)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes analysis
fit2 <- eBayes(fit2)

# Extract results
results_3h <- topTable(
  fit2,
  coef = "H3_vs_Baseline",
  number = Inf,
  adjust.method = "BH"
)

head(results_3h)

# Extract differential expression results

results_3h <- topTable(
  fit2,
  coef = "H3_vs_Baseline",
  number = Inf,
  adjust.method = "BH"
)

results_48h <- topTable(
  fit2,
  coef = "H48_vs_Baseline",
  number = Inf,
  adjust.method = "BH"
)

results_96h <- topTable(
  fit2,
  coef = "H96_vs_Baseline",
  number = Inf,
  adjust.method = "BH"
)

# Count significant probes (adjusted p-value < 0.05)

sum(results_3h$adj.P.Val < 0.05, na.rm = TRUE)

sum(results_48h$adj.P.Val < 0.05, na.rm = TRUE)

sum(results_96h$adj.P.Val < 0.05, na.rm = TRUE)


results_path <- "C:/Users/anopj/OneDrive/Documents/Bioinformatics/results"

dir.create(results_path, showWarnings = FALSE)

write.csv(results_3h,
          file.path(results_path, "DEG_3h_vs_baseline.csv"))

write.csv(results_48h,
          file.path(results_path, "DEG_48h_vs_baseline.csv"))

write.csv(results_96h,
          file.path(results_path, "DEG_96h_vs_baseline.csv"))