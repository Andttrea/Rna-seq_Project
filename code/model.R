# Load libraries
library(recount3)

# Save the variable SRE project data
rse_gene_SRP127181 <- readRDS(file = "processed-data/rse_gene_SRP127181")

# We need to change the data format to use it in the statistical model, so we will cast the data
# First we need to see what type of data we have in each column
colData(rse_gene_SRP127181)[, grepl(
  "^sra_attribute",
  colnames(colData(rse_gene_SRP127181))
)]

# Since we don't have any numeric characters, we will convert the data to factors.


