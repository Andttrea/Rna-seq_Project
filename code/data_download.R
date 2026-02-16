# Load library 
library(recount3)

# Check mouse available projects and save the result in a variable
mouse_project <- available_projects(organism = "mouse") 

# I will use SRP057814 for the analysis because its samples represent a highly controlled experimental design
# Download the data and save it in a variable
proj_info <- subset(mouse_project, project == "SRP057814" &  project_type == "data_sources")

# Create a RangedSummarizedExperiment object with the data
rse_gene_SRP057814 <- create_rse(proj_info) 

