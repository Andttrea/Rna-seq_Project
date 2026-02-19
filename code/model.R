# Load libraries
library(recount3)

# Save the variable SRE project data
rse_gene_SRP057814 <- readRDS(file = "processed-data/rse_gene_SRP057814")

# We need to change the data format to use it in the statistical model, so we will cast the data
# First we need to see what type of data we have in each column
colData(rse_gene_SRP057814)[, grepl(
  "^sra_attribute",
  colnames(colData(rse_gene_SRP057814))
)]

# Since we don't have any numeric characters, we will convert the data to factors.
rse_gene_SRP057814$sra_attribute.cell_type <- factor(tolower(rse_gene_SRP057814$sra_attribute.cell_type))
rse_gene_SRP057814$sra_attribute.chip_antibody <- factor(tolower(rse_gene_SRP057814$sra_attribute.chip_antibody))
rse_gene_SRP057814$sra_attribute.genotype <- factor(tolower(rse_gene_SRP057814$sra_attribute.genotype))
rse_gene_SRP057814$sra_attribute.strain <- factor(rse_gene_SRP057814$sra_attribute.strain)
rse_gene_SRP057814$sra_attribute.source_name <- factor(tolower(rse_gene_SRP057814$sra_attribute.source_name))

