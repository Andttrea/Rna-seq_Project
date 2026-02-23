# Load libraries
library(recount3)
library(edgeR)
library(limma)
library(variancePartition)
library(ggplot2)

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

# Since we need this data in a format that can be used in the statistical model,
# we will convert the data to factors or numeric characters given the type of data
# we have in each column.

# sra_attribute.mutation_status and sra_attribute.source_name are not numerical characters, so we will convert them to factors.
rse_gene_SRP127181$sra_attribute.mutation_status <- factor(
  rse_gene_SRP127181$sra_attribute.mutation_status
)
rse_gene_SRP127181$sra_attribute.source_name <- factor(
  rse_gene_SRP127181$sra_attribute.source_name
)

# The variable sra_attribute.treatment_time has numeric characters, but contains the word "hours" in it. For practical reasons and to be able to
# perform fitExtractVarPartModel later, we will keep it as a factor, just like the previous variables.
rse_gene_SRP127181$sra_attribute.treatment_time <- factor(
  rse_gene_SRP127181$sra_attribute.treatment_time
)

# Checking the data of sra_attribute.treatment, we noticed that there was a typo in the word "ethanol", instead of "ethanol" it was written "ethonol".
# The typo was identified because searching on the internet for "ethonol" does not return any chemical compound.
# Fix the typo and convert the data to factors.
rse_gene_SRP127181$sra_attribute.treatment <- gsub(
  "ethonol",
  "ethanol",
  rse_gene_SRP127181$sra_attribute.treatment
)
rse_gene_SRP127181$sra_attribute.treatment <- factor(
  rse_gene_SRP127181$sra_attribute.treatment
)

# Quality analysis of the data, using the variable assigned_gene_prop
rse_gene_SRP127181$assigned_gene_prop <- rse_gene_SRP127181$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP127181$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP127181$assigned_gene_prop)

# Output of the code above:
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7725  0.7917  0.7984  0.7981  0.8073  0.8140
# Since the quality of the data is good, it is not necessary to remove any sample from the analysis.
# Note: It has been decided that the quality of the data is good because the minimum value of the variable assigned_gene_prop is 0.7725,
# which means that at least 77.25% of the reads were assigned to a CDS.

# Calculate the median levels of gene expression in our samples
gene_means <- rowMeans(assay(rse_gene_SRP127181, "counts"))
summary(gene_means)
# Output of the code above:
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.000e+00 0.000e+00 1.850e-01 5.463e+02 4.520e+01 2.655e+05
# It is noticed that the mean of gene expression is higher than the median, which means that there are some genes with very high levels of expression that are affecting the mean.
# Due to this, it is necessary to filter the data to remove genes with very low levels of expression, because they are not informative for the analysis and can affect the results.

# Filter the data
# Note: We will use the function filterByExpr from the edgeR package
# to filter the data, which uses a method based on counts per million (CPM) to filter the data.

# Create a DGEList object with the data
dge <- DGEList(counts = assay(rse_gene_SRP127181, "counts"))

# Create a vector with the group information of the samples
# Note: The group information is created by concatenating the variables sra_attribute.mutation_status, sra_attribute.treatment and sra_attribute.treatment_time, because these variables are the ones that we will use in the statistical model.
group <- with(
  colData(rse_gene_SRP127181),
  paste0(
    sra_attribute.mutation_status,
    "_",
    sra_attribute.treatment,
    "_",
    sra_attribute.treatment_time
  )
)

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
# This means that we have removed 66.69% of the genes from the original data, which is a good filtering percentage, because we have removed the genes with very low levels of expression that are not informative for the analysis and can affect the results.

# Normalizing the data
# Note: Normalization is necessary to remove technical variation in the data and to make the data comparable between samples.

# Note: The ColData is being used to map samples to metadata to link gene expression with experimental conditions for accurate statistical analysis.
dge <- DGEList(
  counts = assay(rse_gene_SRP127181_filtered, "counts"),
  genes = rowData(rse_gene_SRP127181_filtered),
  samples = colData(rse_gene_SRP127181_filtered)
)
dge <- calcNormFactors(dge)

# Exploring the data with the variancePartition package to see how the different variables are affecting the variance of the data.
# Transform the data to log2 counts per million (logCPM) to establish variance
# Note: VariancePartition requires the data to be in log2 counts per million (logCPM) format
exp_log2 <- voom(dge)

# Define the variables that we want to evaluate
# Note: Using the variable (1/variable) to evaluate the variance of each variable independently
form_var <- ~ (1 | sra_attribute.mutation_status) +
  (1 | sra_attribute.treatment) +
  (1 | (sra_attribute.treatment_time))

# Calculate the partition of the variance
# vGene$E give the logCPM values of the data, form give the variables that we want to evaluate
var_Part <- fitExtractVarPartModel(
  exp_log2$E,
  form_var,
  colData(rse_gene_SRP127181_filtered)
)

# Plot the partition of the variance to see how the different variables are affecting the variance of the data.
plotVarPart(sortCols(var_Part))
# As we can see in the plot, the variable sra_attribute.mutation_status is the one that is affecting the variance of the data the most,
# affecting around 50% of the variance, the other variables do not affect the variance of the data as much.
# To know more about the variance between the different variables, a PCA is being used

# Extract normalized counts using edgeR's cpm function
counts_norm <- edgeR::cpm(dge, log = TRUE)

# Perform PCA
pca <- prcomp(t(counts_norm), scale. = TRUE)

# Create data frame for plotting
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Mutation = rse_gene_SRP127181_filtered$sra_attribute.mutation_status
)

# Plot with color by mutation only
ggplot(pca_df, aes(x = PC1, y = PC2, color = Mutation)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCA by Mutation Status",
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 2), "%)")
  )

# It is noticed that the different mutation statuses are separated in the PCA plot, only the Wild Type is divided into two groups,
# this could be due to another variable affecting the variance of the data or to batch effect.
# To design our model matrix, we are considering not only the variable sra_attribute.mutation_status, but also the variable sra_attribute.treatment_time because it is the second variable that is affecting the variance.
# This is with the purpose of having a good model.

# Before creating the model we need to reorder the levels

# Note: If we observe the levels of our variable sra_attribute.mutation_status, we see the following:
# levels(rse_gene_SRP127181$sra_attribute.mutation_status)
#[1] "ESR1 D538G" "ESR1 WT"    "ESR1 Y537S"
# If we create the model with that order in the levels, then we would have our mutation "ESR1 D538G" as the base, which
# would be what we compare with the mutation "ESR1 Y537S" and the WT. In this case, it would be more interesting to observe
# the difference between the WT and the two mutations, so we would need "ESR1 WT" to be our base.

levels(rse_gene_SRP127181_filtered$sra_attribute.mutation_status)
rse_gene_SRP127181_filtered$sra_attribute.mutation_status <- relevel(
  rse_gene_SRP127181_filtered$sra_attribute.mutation_status,
  "ESR1 WT",
  "ESR1 D538G",
  "ESR1 Y537S"
)
# The order of our levels is now:
# [1] "ESR1 WT"    "ESR1 D538G" "ESR1 Y537S"

# Create model
mod <- model.matrix(
  ~ sra_attribute.mutation_status +
    sra_attribute.treatment_time +
    assigned_gene_prop,
  data = colData(rse_gene_SRP127181_filtered)
)
colnames(mod)
# Output of the code above:
#[1] "(Intercept)"                             "sra_attribute.mutation_statusESR1 D538G" "sra_attribute.mutation_statusESR1 Y537S"
#[4] "sra_attribute.treatment_time4"           "assigned_gene_prop"

# Note: Thanks to "sra_attribute.treatment_time4" we can realize that the time used as reference is 24 hours, in this case
# we will not reorder it.

# Save variables in processed-data folder for later use

saveRDS(
  rse_gene_SRP127181_filtered,
  file = "processed-data/rse_gene_SRP127181_filtered"
)
saveRDS(dge, file = "processed-data/dge")
saveRDS(mod, file = "processed-data/model_matrix")
