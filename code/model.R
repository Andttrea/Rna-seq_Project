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

# Output of the code above:
# DataFrame with 27 rows and 4 columns
#           sra_attribute.mutation_status sra_attribute.source_name
#                             <character>               <character>
#SRR6392012                       ESR1 WT                      T47D
#SRR6392013                       ESR1 WT                      T47D
#SRR6392014                       ESR1 WT                      T47D
#SRR6392015                       ESR1 WT                      T47D
#SRR6392016                       ESR1 WT                      T47D
#...                                  ...                       ...
#SRR6392034                    ESR1 D538G                       TDG
#SRR6392035                    ESR1 D538G                       TDG
#SRR6392036                    ESR1 Y537S                       TYS
#SRR6392037                    ESR1 Y537S                       TYS
#SRR6392038                    ESR1 Y537S                       TYS
#           sra_attribute.treatment_time sra_attribute.treatment
#                            <character>             <character>
#SRR6392012                      4 hours                 ethonol
#SRR6392013                      4 hours                 ethonol
#SRR6392014                      4 hours                 ethonol
#SRR6392015                      4 hours           10nM estrogen
#SRR6392016                      4 hours           10nM estrogen
#...                                 ...                     ...
#SRR6392034                     24 hours           10nM estrogen
#SRR6392035                     24 hours           10nM estrogen
#SRR6392036                     24 hours           10nM estrogen
#SRR6392037                     24 hours           10nM estrogen
#SRR6392038                     24 hours           10nM estrogen

# Due we need this data in a format that can be used in the statistical model,
#  we will convert the data to factors or numeric characters given the type of data
#  we have in each column.

# sra_attribute.mutation_status and sra_attribute.source_name aren´t numerical characters , so it will convert to factors.
rse_gene_SRP127181$sra_attribute.mutation_status <- factor(rse_gene_SRP127181$sra_attribute.mutation_status)
rse_gene_SRP127181$sra_attribute.source_name <- factor(rse_gene_SRP127181$sra_attribute.source_name)

# The variable sra_attribute.treatment_time has numeric characters, but have the word "hours" in it, so we will remove the word "hours" and then convert the data to numeric characters.
rse_gene_SRP127181$sra_attribute.treatment_time <- gsub(" hours", "", rse_gene_SRP127181$sra_attribute.treatment_time)
rse_gene_SRP127181$sra_attribute.treatment_time <- as.numeric(rse_gene_SRP127181$sra_attribute.treatment_time)

# Checking the data of sra_attribute.treatment, I noticed that there is a typo in the word "ethanol", instead of "ethanol" it is written "ethonol".
# I relize that "ethonol" was a typo because searching in internet for "ethonol" does not give any quemical compound. 
# I decide to correct the typo and then convert the data to a factor like the other variables.
rse_gene_SRP127181$sra_attribute.treatment <- gsub("ethonol", "ethanol", rse_gene_SRP127181$sra_attribute.treatment)
rse_gene_SRP127181$sra_attribute.treatment <- factor(rse_gene_SRP127181$sra_attribute.treatment)
 

