# Load libraries
library(limma)

# Guardando anteriores variables 
rse_gene_SRP127181_filtered <- readRDS(file = "processed-data/rse_gene_SRP127181_filtered")
dge <- readRDS(file = "processed-data/dge")
model <- readRDS(file = "processed-data/model_matrix")

# Voom transformation
vGene <- voom(dge, model, plot = TRUE)

# Fit linear model
eb_results <- eBayes(lmFit(vGene, model))
# Get the top differentially expressed genes
DE_results <- topTable(eb_results, coef = 2, number = nrow(rse_gene_SRP127181_filtered), sort.by = "none")
dim(DE_results)



head(DE_results)
table(DE_results$adj.P.Val < 0.05)
plotMA(eb_results, coef = 2)

volcanoplot(eb_results, coef = 2, highlight = 4, names = DE_results$gene_name)
