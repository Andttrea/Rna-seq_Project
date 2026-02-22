# Load libraries
library(recount3)
library(edgeR)

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

# Checking the data of sra_attribute.treatment, it was notice that there was a typo in the word "ethanol", instead of "ethanol" it is written "ethonol".
# The typo was known because searching in internet for "ethonol" does not give any quemical compound. 
# Fix the typo and convert the data to factors.
rse_gene_SRP127181$sra_attribute.treatment <- gsub("ethonol", "ethanol", rse_gene_SRP127181$sra_attribute.treatment)
rse_gene_SRP127181$sra_attribute.treatment <- factor(rse_gene_SRP127181$sra_attribute.treatment)
 
# Quality analysis of the data, using the variable assigned_gene_prop
rse_gene_SRP127181$assigned_gene_prop <- rse_gene_SRP127181$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP127181$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP127181$assigned_gene_prop)

# Outut of the code above:
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7725  0.7917  0.7984  0.7981  0.8073  0.8140
# Due the quality of the data is good, is not necessary to remove any sample from the analysis. 
# Note: It is been decided that the quality of the data is goood because the minimum value of the variable assigned_gene_prop is 0.7725, 
# which means that at least 77.25% of the reads were assigned to a CDS.

# Calculate the median levels of gene expression in our samples 
gene_means <- rowMeans(assay(rse_gene_SRP127181, "counts"))
summary(gene_means)
# Output of the code above:
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 0.000e+00 1.850e-01 5.463e+02 4.520e+01 2.655e+05 
# It is notice that the mean of gene expression is higher than the median, which means that there are some genes with very high levels of expression that are affecting the mean.
# Due to this, it is necessary to filter the data to remove the genes with very low levels of expression, because they are not informative for the analysis and can affect the results.

# Filter the data
# Note:We will use the function filterByExpr from the edgeR package 
# to filter the data, which uses a method based on the counts per million (CPM) to filter the data.

# Create a DGEList object with the data
dge <- DGEList(counts = assay(rse_gene_SRP127181, "counts"))

# Create a vector with the group information of the samples
# Note: The group information is created by concatenating the variables sra_attribute.mutation_status, sra_attribute.treatment and sra_attribute.treatment_time, because these variables are the ones that we will use in the statistical model.
group <- with(colData(rse_gene_SRP127181), paste0(sra_attribute.mutation_status, "_", sra_attribute.treatment, "_", sra_attribute.treatment_time))

# Filter the data using the function filterByExpr
keep <- filterByExpr(dge, group = group)

# Filter the original RangedSummarizedExperiment 
rse_gene_SRP127181_filtered <- rse_gene_SRP127181[keep, ]

# Check the dimensions of the original and filtered data
dim(rse_gene_SRP127181)
# Output of the code above:
# [1] 63856    27
dim(rse_gene_SRP127181_filtered)
# Output of the code above:
# [1] 21272    27

# Check the ratio of the number of genes in the original and filtered data
round(nrow(rse_gene_SRP127181_filtered) / nrow(rse_gene_SRP127181) * 100, 2)
# Output of the code above:
# [1] 33.31



