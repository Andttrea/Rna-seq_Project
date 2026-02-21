# Load library 
library(recount3)

# Check human available projects and save the result in a variable
human_project <- available_projects() 

# I will use SRP127181 for the analysis because its samples represent a highly controlled experimental design
# Also have various SRA attributes that we can use in the statistical model. 

# Download the data and save it in a variable
proj_info <- subset(human_project, project == "SRP127181" &  project_type == "data_sources")

# Create a RangedSummarizedExperiment object with the data
rse_gene_SRP127181 <- create_rse(proj_info) 
rse_gene_SRP127181


# Compute read counts and save the result in the assay slot of the RangedSummarizedExperiment object
assay(rse_gene_SRP127181, "counts") <- compute_read_counts(rse_gene_SRP127181)

# Expand the SRA attributes and save the result in the colData slot of the RangedSummarizedExperiment object
rse_gene_SRP127181 <- expand_sra_attributes(rse_gene_SRP127181) 


# Save the RangedSummarizedExperiment object in an RDS file
saveRDS(create_rse(proj_info), file = "raw-data/raw_rse_gene_SRP127181") # Raw RSE
saveRDS(rse_gene_SRP127181, file = "processed-data/rse_gene_SRP127181") # Processed RSE
