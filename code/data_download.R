# Load library 
library(recount3)

# Check mouse available projects and save the result in a variable
mouse_project <- available_projects(organism = "mouse") 

# I will use SRP057814 for the analysis because its samples represent a highly controlled experimental design
# Also 

# Download the data and save it in a variable
proj_info <- subset(mouse_project, project == "SRP057814" &  project_type == "data_sources")

# Create a RangedSummarizedExperiment object with the data
rse_gene_SRP057814 <- create_rse(proj_info) 
rse_gene_SRP057814


# Compute read counts and save the result in the assay slot of the RangedSummarizedExperiment object
assay(rse_gene_SRP057814, "counts") <- compute_read_counts(rse_gene_SRP057814)

# Expand the SRA attributes and save the result in the colData slot of the RangedSummarizedExperiment object
rse_gene_SRP057814 <- expand_sra_attributes(rse_gene_SRP057814) 