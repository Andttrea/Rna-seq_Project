# Finally, we will examine the treatment time variable to see if there are more genes expressed over time.
# We perform the same process as we did for the mutationsGrison between Wild Type and Mutant D538G
model <- readRDS(file = "processed-data/model_matrix")

# Voom transformation
vGene <- voom(dge, model, plot = TRUE)

# Fit linear model
eb_results <- eBayes(lmFit(vGene, model))
# Get the top differentially expressed genes para la comparación entre el Wild Type y el Mutant D538G
DE_results_D538G <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP127181_filtered),
  sort.by = "none"
)
# We get the table dimensions to check if we didn't lose any genes in the analysis process
dim(DE_results_D538G)
# Output of the code above:
# [1] 21272    27
# Observamos que no perdimos ningún gen

# Observamos DE_results
head(DE_results_D538G)
# Output of the code above:
#                   source type bp_length phase           gene_id              gene_type     gene_name level
#ENSG00000227232.5 HAVANA gene      1351    NA ENSG00000227232.5 unprocessed_pseudogene        WASH7P     2
#ENSG00000238009.6 HAVANA gene      3726    NA ENSG00000238009.6                lincRNA  RP11-34P13.7     2
#ENSG00000233750.3 HAVANA gene      3812    NA ENSG00000233750.3   processed_pseudogene        CICP27     1
#ENSG00000268903.1 HAVANA gene       755    NA ENSG00000268903.1   processed_pseudogene RP11-34P13.15     2
#ENSG00000269981.1 HAVANA gene       284    NA ENSG00000269981.1   processed_pseudogene RP11-34P13.16     2
#ENSG00000241860.6 HAVANA gene      6195    NA ENSG00000241860.6   processed_transcript RP11-34P13.13     2
#                           havana_gene               tag       logFC   AveExpr          t      P.Value
#ENSG00000227232.5 OTTHUMG00000000958.1              <NA> -0.15453287  1.882552 -1.3539947 1.869565e-01
#ENSG00000238009.6 OTTHUMG00000001096.2 overlapping_locus -0.04509432 -1.343374 -0.1727149 8.641630e-01
#ENSG00000233750.3 OTTHUMG00000001257.3    pseudo_consens -1.04283697 -2.372336 -3.6438460 1.125553e-03
#ENSG00000268903.1 OTTHUMG00000182518.2              <NA> -1.40238502 -1.288555 -5.8482074 3.151641e-06
#ENSG00000269981.1 OTTHUMG00000182738.2              <NA> -0.79611373 -1.671208 -2.9636278 6.277981e-03
#ENSG00000241860.6 OTTHUMG00000002480.3        ncRNA_host -0.08292388 -2.005890 -0.3208142 7.508197e-01
#                     adj.P.Val          B
#ENSG00000227232.5 2.430446e-01 -6.4310394
#ENSG00000238009.6 8.912715e-01 -6.5727263
#ENSG00000233750.3 2.622715e-03 -0.9341409
#ENSG00000268903.1 1.513357e-05  4.5514185
#ENSG00000269981.1 1.207461e-02 -2.6452264
#ENSG00000241860.6 7.948361e-01 -6.3517654

# Filter genes with P-value less than 0.05
table(DE_results_D538G$adj.P.Val < 0.05)
# Output of the code above:
#FALSE  TRUE
# 8123 13149
# Se observa que hay 13149 genes que son significativamente diferentes entre el Wild Type y el Mutant D538G.

# Create an MA plot
plotMA(eb_results, coef = 2)

# Create a volcano plot, highlighting the 4 most significant genes
volcanoplot(
  eb_results,
  coef = 2,
  highlight = 4,
  names = DE_results_D538G$gene_name
)

# Lo anterior solo fue hecho para la comparación entre el Wild Type y el Mutant D538G
# Necesitaremos comparar la otra mutación Y537S con el Wild Type, para eso haremos lo mismo que hicimos para D538G

DE_results_Y537S <- topTable(
  eb_results,
  coef = 3,
  number = nrow(rse_gene_SRP127181_filtered),
  sort.by = "none"
)
dim(DE_results_Y537S)
# Output of the code above:
#[1] 21272    16
# No genes were lost in the analysis process

head(DE_results_Y537S)
# Output of the code above:
#   source type bp_length phase           gene_id              gene_type     gene_name level
#ENSG00000227232.5 HAVANA gene      1351    NA ENSG00000227232.5 unprocessed_pseudogene        WASH7P     2
#ENSG00000238009.6 HAVANA gene      3726    NA ENSG00000238009.6                lincRNA  RP11-34P13.7     2
#ENSG00000233750.3 HAVANA gene      3812    NA ENSG00000233750.3   processed_pseudogene        CICP27     1
#ENSG00000268903.1 HAVANA gene       755    NA ENSG00000268903.1   processed_pseudogene RP11-34P13.15     2
#ENSG00000269981.1 HAVANA gene       284    NA ENSG00000269981.1   processed_pseudogene RP11-34P13.16     2
#ENSG00000241860.6 HAVANA gene      6195    NA ENSG00000241860.6   processed_transcript RP11-34P13.13     2
#                           havana_gene               tag      logFC   AveExpr         t      P.Value    adj.P.Val          B
#ENSG00000227232.5 OTTHUMG00000000958.1              <NA> -0.4608210  1.882552 -3.383614 2.199146e-03 0.0043286966 -2.5505247
#ENSG00000238009.6 OTTHUMG00000001096.2 overlapping_locus -0.4774108 -1.343374 -1.477651 1.510652e-01 0.1975689262 -5.6134246
#ENSG00000233750.3 OTTHUMG00000001257.3    pseudo_consens -1.4389688 -2.372336 -4.176659 2.766436e-04 0.0006618785  0.3156886
#ENSG00000268903.1 OTTHUMG00000182518.2              <NA> -1.2966435 -1.288555 -4.677415 7.237972e-05 0.0001974684  1.4198915
#ENSG00000269981.1 OTTHUMG00000182738.2              <NA> -0.9573642 -1.671208 -3.092189 4.575970e-03 0.0084284378 -2.4875362
#ENSG00000241860.6 OTTHUMG00000002480.3        ncRNA_host -0.4846160 -2.005890 -1.541983 1.347107e-01 0.1782067787 -5.3617502

table(DE_results_Y537S$adj.P.Val < 0.05)
# Output of the code above:
# FALSE  TRUE
# 7358 13914

# Create an MA plot
plotMA(eb_results, coef = 3)

# Create a volcano plot, highlighting the 4 most significant genes
volcanoplot(
  eb_results,
  coef = 3,
  highlight = 4,
  names = DE_results_Y537S$gene_name
)

# Por ultimo, observaremos la variable del tiempo de tratamiento, para ver si hay más genes expresados a través del tiempo.
# Realizamos el mismo proceso que hicimos para las mutaciones

DE_results_time <- topTable(
  eb_results,
  coef = 4,
  number = nrow(rse_gene_SRP127181_filtered),
  sort.by = "none"
)
dim(DE_results_time)
# Output of the code above:
# [1] 21272    16

head(DE_results_time)
# Output of the code above:
#    source type bp_length phase           gene_id              gene_type     gene_name level
#ENSG00000227232.5 HAVANA gene      1351    NA ENSG00000227232.5 unprocessed_pseudogene        WASH7P     2
#ENSG00000238009.6 HAVANA gene      3726    NA ENSG00000238009.6                lincRNA  RP11-34P13.7     2
#ENSG00000233750.3 HAVANA gene      3812    NA ENSG00000233750.3   processed_pseudogene        CICP27     1
#ENSG00000268903.1 HAVANA gene       755    NA ENSG00000268903.1   processed_pseudogene RP11-34P13.15     2
#ENSG00000269981.1 HAVANA gene       284    NA ENSG00000269981.1   processed_pseudogene RP11-34P13.16     2
#ENSG00000241860.6 HAVANA gene      6195    NA ENSG00000241860.6   processed_transcript RP11-34P13.13     2
#                           havana_gene               tag       logFC   AveExpr          t     P.Value  adj.P.Val         B
#ENSG00000227232.5 OTTHUMG00000000958.1              <NA> -0.03177697  1.882552 -0.3887589 0.700501487 0.80553672 -6.775980
#ENSG00000238009.6 OTTHUMG00000001096.2 overlapping_locus  0.29103545 -1.343374  1.4753110 0.151688896 0.29043440 -5.030903
#ENSG00000233750.3 OTTHUMG00000001257.3    pseudo_consens  0.35182339 -2.372336  1.7075675 0.099185176 0.21389334 -4.528882
#ENSG00000268903.1 OTTHUMG00000182518.2              <NA>  0.47367915 -1.288555  2.7965083 0.009399936 0.03916088 -2.738240
#ENSG00000269981.1 OTTHUMG00000182738.2              <NA>  0.60677409 -1.671208  3.0129370 0.005563996 0.02679586 -2.207284
#ENSG00000241860.6 OTTHUMG00000002480.3        ncRNA_host  0.24237037 -2.005890  1.2641803 0.216963649 0.37236367 -5.152552

table(DE_results_time$adj.P.Val < 0.05)
# Output of the code above:
#FALSE  TRUE
#15658  5614

# Plots
plotMA(eb_results, coef = 4)
volcanoplot(
  eb_results,
  coef = 4,
  highlight = 4,
  names = DE_results_time$gene_name
)

# Finally, we will create a heatmap with the 50 most significant genes

# Select the top 50 genes with the lowest overall p-value (F-statistic)
# eb_results$F.p.value tests if ANY coefficient in the model is non-zero
top_indices <- order(eb_results$F.p.value)[1:50]
global_heatmap_matrix <- vGene$E[top_indices, ]

# 2. Assign human-readable gene names to the rows
# Replacing Ensembl IDs with Gene Symbols for better interpretability
rownames(global_heatmap_matrix) <- rse_gene_SRP127181_filtered$gene_name[
  top_indices
]

# 3. Prepare sample metadata for the heatmap annotation
# We include all variables to see how they influence the clustering
metadata_annotation <- as.data.frame(colData(rse_gene_SRP127181_filtered)[,
  c(
    "sra_attribute.mutation_status",
    "sra_attribute.treatment_time",
    "sra_attribute.treatment"
  )
])

# Rename columns for a cleaner legend in the plot
colnames(metadata_annotation) <- c("Mutation", "Time", "Treatment")

pheatmap(
  global_heatmap_matrix,
  cluster_rows = TRUE, # Groups genes with similar expression trends
  cluster_cols = TRUE, # Groups samples with similar transcriptional profiles
  show_colnames = FALSE, # Hide individual sample IDs for a cleaner look
  annotation_col = metadata_annotation,
  scale = "row", # Standardizes expression per gene (Z-score)
  main = "Top 50 DE Genes: Global Experiment View",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
