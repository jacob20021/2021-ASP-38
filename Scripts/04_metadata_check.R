# Check sample time points
table(meta_pmn$`time point:ch1`)

# Check participant IDs
table(meta_pmn$`subject id:ch1`)

# Check number of samples per time point
table(
  meta_pmn$`subject id:ch1`,
  meta_pmn$`time point:ch1`
)

# Check metadata structure
str(meta_pmn)