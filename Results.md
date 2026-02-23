# Diferential Expression Analysis

## Experimental Variables and Background

### Mutation Status (mutant_status)

The variable `mutant_status` describes the genotype of the ESR1 (Estrogen Receptor Alpha)
gene in the T47D breast cancer cell line used in this study. Three genotypes were analyzed:

**ESR1 Wild Type (WT):**
    The normal, unmutated form of the estrogen receptor alpha gene. In its
    wild-type state, ESR1 encodes a ligand-dependent transcription factor
    that requires estrogen binding to activate target gene expression.
    WT ESR1 is the standard against which mutation effects are measured
    in breast cancer research.

**ESR1 Y537S:**
    A somatic point mutation in the ligand-binding domain of ESR1, where
    tyrosine (Y) at position 537 is replaced by serine (S). This mutation
    stabilizes the active conformation of the receptor, resulting in
    constitutive, estrogen-independent transcriptional activity. Y537S is
    one of the two most frequent ESR1 mutations found in metastatic
    ER-positive breast cancer and is associated with resistance to
    endocrine therapies such as aromatase inhibitors.

**ESR1 D538G:**
    A somatic point mutation in the ligand-binding domain of ESR1, where
    aspartic acid (D) at position 538 is replaced by glycine (G). Similar
    to Y537S, this mutation confers ligand-independent activation of the
    estrogen receptor, enabling constitutive gene expression without
    estrogen. D538G is frequently detected in metastatic hormone
    receptor-positive breast cancers that have acquired resistance to
    standard endocrine treatments.

### Treatment Conditions

The experimental treatments consisted of ethanol as the vehicle control
and 10 nM estrogen (E2) as the active treatment. Ethanol serves as
the solvent in which estrogen is dissolved, so it is used as a baseline
to distinguish the biological effects of estrogen from any solvent-related
artifacts. Samples were collected at two time points: 4 hours and 24 hours
after treatment, allowing the observation of both early and sustained
transcriptional responses to estrogen across the different ESR1 genotypes.

### Rationale

Comparing wild-type and mutant ESR1 cells under both vehicle and estrogen
conditions at different time points provides a framework to understand how
these clinically relevant mutations reprogram the transcriptome. Since Y537S
and D538G mutations confer estrogen-independent receptor activity, this design
allows the identification of genes that are aberrantly activated by the mutant
receptor even in the absence of hormone, and reveals how the transcriptional
response to estrogen differs between normal and mutant receptor contexts.

---

## Dataset Description

The data used in this analysis was obtained from the Sequence Read Archive
project SRP127181 (GEO accession GSE108304). This dataset contains RNA-seq
data from 27 human samples corresponding to the T47D breast cancer cell line
and its genome-edited variants harboring ESR1 Y537S (TYS) and ESR1 D538G (TDG)
mutations. Hormone-depleted cells were treated with either 10 nM estrogen or
ethanol (vehicle control) for 4 or 24 hours, with 3 biological replicates per
condition. Sequencing was performed on the Illumina HiSeq 4000 platform. After
quality control and filtering of lowly expressed genes, 21,272 genes were retained
for differential expression analysis using the limma-voom pipeline.

---

## Results

### 1. Mean-Variance Trend (Voom Transformation)

![Mean-Variance Trend](plots/Mean_variance.svg)

The mean-variance plot displays the relationship between the average expression
level of each gene (x-axis, in log2 counts per million) and its standard deviation
(y-axis). This plot is generated during the voom transformation, which models
the mean-variance relationship to assign precision weights to each observation,
enabling the use of linear models on RNA-seq count data. The red curve represents
the fitted trend that captures how variance changes with expression level.

In our data, the trend line shows a smooth decreasing pattern from left to right,
indicating that lowly expressed genes exhibit higher variability while highly
expressed genes show more stable variance. The points are evenly distributed around
the fitted trend without major deviations, which confirms that the voom transformation
successfully stabilized the variance across the full range of gene expression. This
result indicates that the data is well-suited for downstream linear modeling and that
no major technical artifacts or batch effects are distorting the analysis. The broad
dynamic range observed (approximately -5 to 15 log2 CPM) reflects the transcriptomic
diversity captured across the 27 samples.

---

### 2. MA Plot: WT vs D538G

![MA Plot WT vs D538G](plots/plotMA_WT_D538G.svg)

An MA plot visualizes differential expression by plotting the average expression of
each gene (A, x-axis) against its log-fold change (M, y-axis) between two conditions.
Genes centered around a log-fold change of zero show no difference between conditions,
while genes displaced above or below zero are upregulated or downregulated, respectively.
Red points indicate statistically significant genes (adjusted p-value < 0.05).

In the comparison between ESR1 Wild Type and D538G mutant cells, the majority of genes
cluster around a logFC of 0, as expected. However, a substantial number of genes show
positive or negative fold changes, indicating widespread transcriptional reprogramming
caused by the D538G mutation. A total of 13,149 genes were found to be differentially
expressed (FDR < 0.05), representing approximately 61.8% of the analyzed transcriptome.
This suggests that the D538G mutation induces broad alterations in gene expression,
consistent with the constitutive activation of the estrogen receptor driving large-scale
changes in downstream signaling pathways.

---

### 3. MA Plot: WT vs Y537S

![MA Plot WT vs Y537S](plots/plotMA_WT_Y537S.svg)

In the comparison between ESR1 Wild Type and Y537S mutant cells, the data follows a
similar overall structure to the D538G comparison, with most genes centered at logFC
of 0. Nevertheless, a large proportion of genes deviate from zero, indicating significant
differential expression. The Y537S comparison identified 13,914 differentially expressed
genes (FDR < 0.05), which is slightly more than the D538G comparison. This is consistent
with published reports showing that Y537S exhibits a more pronounced and unique
transcriptional phenotype compared to D538G. The spread of fold changes in both positive
and negative directions indicates that the Y537S mutation both activates and represses
distinct sets of genes relative to the wild-type receptor.

---

### 4. MA Plot: Treatment Time (4h vs 24h)

![MA Plot Treatment Time](plots/plotMA_time.svg)

This MA plot compares gene expression between the 4-hour and 24-hour treatment time
points across all genotypes. The majority of genes are centered near logFC of 0,
indicating that many genes maintain stable expression regardless of treatment duration.
However, 5,614 genes were identified as differentially expressed between the two time
points (FDR < 0.05). This is considerably fewer than the mutation comparisons, suggesting
that while time influences transcriptional dynamics, the ESR1 mutation status has a
stronger overall effect on gene expression. The genes showing time-dependent changes
likely represent secondary estrogen response genes that require prolonged signaling
to become activated or repressed.

---

### 5. Volcano Plot: WT vs D538G

![Volcano Plot WT vs D538G](plots/volcano_WT_D538G.svg)

A volcano plot combines statistical significance (y-axis, -log10 p-value) with biological
effect size (x-axis, log2 fold change) to identify genes that are both statistically
significant and biologically meaningful. Genes in the upper corners of the plot represent
the most relevant candidates, as they have both large fold changes and high statistical
confidence. The genes highlighted in blue represent the four most statistically significant
genes in this comparison.

In the WT vs D538G volcano plot, the distribution of points reveals a large number of genes with significant p-values on both sides of the fold-change axis, confirming the widespread transcriptional impact of the D538G mutation. Four genes are highlighted as the most statistically significant hits: PIP, UPK1A, and SDC2 are located on the left side of the plot, indicating they are significantly downregulated in D538G mutant cells, while MAGEB2 appears on the right side, indicating it is significantly upregulated.

**PIP (Prolactin-Induced Protein)** is a well-established luminal breast cancer marker and estrogen-responsive gene; its downregulation may reflect a shift away from a differentiated luminal phenotype driven by the constitutive activation of the mutant receptor. **UPK1A (Uroplakin 1A)** is associated with epithelial integrity and its downregulation suggests alterations in epithelial identity and membrane-associated processes induced by the mutation. **SDC2 (Syndecan-2)**, a cell-surface heparan sulfate proteoglycan involved in cell–matrix interactions, adhesion, and growth factor signaling, is also suppressed, which may reflect a reorganization of the tumor microenvironment interactions in D538G mutant cells. In contrast, **MAGEB2 (Melanoma Antigen Family B2)**, a cancer-testis antigen, is significantly upregulated, consistent with its known roles in promoting oncogenic phenotypes including enhanced proliferation, immune evasion, and tumor aggressiveness.

Together, the downregulation of PIP, UPK1A, and SDC2 alongside the upregulation of MAGEB2 supports a model in which the D538G mutation drives transcriptional reprogramming that promotes a more aggressive, ligand-independent tumor phenotype.

---

### 6. Volcano Plot: WT vs Y537S

![Volcano Plot WT vs Y537S](plots/volcano_WT_Y537S.svg)

The volcano plot for the **WT vs Y537S** comparison shows a broad distribution of differentially expressed genes, with many reaching very high statistical significance (−log₁₀ p-value > 20). Similar to the D538G comparison, genes are displaced in both directions along the fold change axis, but the Y537S mutation appears to produce an even more pronounced transcriptional shift. The four genes highlighted in blue represent the most statistically significant changes in this dataset: **UPK1A**, **NANOS1**, and **DNAJA4**, which are downregulated, and **BCHE**, which is markedly upregulated.

Among these, **UPK1A (Uroplakin 1A)** shows the strongest negative fold change, indicating substantial repression in **Y537S** mutant cells. Although classically associated with epithelial differentiation and membrane structure, its downregulation here may reflect alterations in epithelial identity and cell surface organization driven by constitutive ER signaling. **NANOS1**, an RNA-binding protein involved in post-transcriptional regulation and cell fate control, is also significantly suppressed, suggesting potential rewiring of RNA regulatory networks in the mutant context. **DNAJA4**, a member of the Hsp40 (DnaJ) chaperone family, participates in protein folding and stress response pathways; its downregulation may indicate altered proteostasis or stress adaptation mechanisms in Y537S cells.

In contrast, **BCHE (Butyrylcholinesterase)** is highly upregulated and displays one of the strongest positive fold changes in the plot. Although not traditionally considered a canonical estrogen-responsive gene, BCHE overexpression has been reported in several malignancies and may reflect metabolic or microenvironmental adaptations associated with aggressive tumor phenotypes.

---

### 7. Volcano Plot: Treatment Time (4h vs 24h)

![Volcano Plot Treatment Time](plots/volcano_time.svg)

The volcano plot for the time comparison (4 hours vs. 24 hours as base) shows a more compact distribution compared to the mutation comparisons, with fewer genes reaching extreme statistical significance. This is consistent with the lower number of differentially expressed genes identified in this comparison (5,614). The four genes highlighted in blue represent the strongest time-dependent changes across all conditions. In this plot, where 24 hours serves as the reference group, **UPK1A** and **CGA** are located on the right side (positive fold change), indicating they are more highly expressed at 4 hours and subsequently downregulated as the treatment progresses to 24 hours. Conversely, **FGG** and **ANXA9** are on the left side (negative fold change), showing they are significantly upregulated at the 24-hour mark compared to the early 4-hour response.

The fact that **UPK1A** (Uroplakin 1A) and **CGA** (Glycoprotein Hormones Alpha Polypeptide) appear upregulated at 4h relative to 24h suggests they are "early-peak" genes. Their expression is rapidly induced by the initial stimulus but is not sustained, potentially due to negative feedback loops or the transient nature of early estrogen receptor (ER) recruitment to their promoters. On the other hand, the accumulation of **FGG (Fibrinogen Gamma Chain)** and **ANXA9 (Annexin A9)** at 24 hours points toward a late-stage transcriptional program. ANXA9, involved in calcium-dependent membrane signaling, may represent a secondary response that requires prolonged ER activity or the synthesis of intermediate transcription factors to reach its peak expression.

These time-responsive genes likely include estrogen-regulated genes that show a delayed transcriptional response, requiring sustained receptor signaling over 24 hours to reach full activation or repression. This temporal pattern is expected for secondary target genes that depend on intermediate transcription factors or chromatin remodeling events initiated by the early estrogen response. The dynamic behavior of UPK1A, which is suppressed by ESR1 mutations but shows a strong early-induction peak at 4 hours, highlights the complex temporal regulation of the ER-driven transcriptome and how the duration of signaling fundamentally shifts the cellular identity.

---

### 8. Heatmap: Top 50 Differentially Expressed Genes

![Heatmap](plots/heatmap_mutation_time.svg)

A heatmap represents gene expression values as a color gradient, where each row
corresponds to a gene and each column to a sample. Colors indicate the relative
expression level: red denotes upregulation and blue denotes downregulation, with values
standardized per gene using Z-scores. Hierarchical clustering of both rows and columns
groups genes and samples with similar expression patterns, revealing underlying
biological structure in the data.

In our heatmap reveals a clear hierarchy in mutational potency, where the **Y537S** variant (green label) exerts the most profound and dominant impact on transcriptional reprogramming, displaying significantly more intense and defined red and blue color contrasts than the other conditions. For its part, the **D538G** mutation (pink label) also achieves constitutive activation of the estrogen receptor independent of the hormone, but with a slightly more moderate magnitude of change and expression profiles that, while similar, show distinct nuances in the saturation of the gene clusters. While **Y537S** stands out as the biochemically most potent and aggressive mutation by more effectively stabilizing the active conformation of the receptor, D538G represents an alternative pathway of resistance that, although visually less extreme, shares the ability to decouple tumor growth from estrogen, consolidating both as the primary drivers of the persistent oncogenic signaling observed across the expression blocks in the plot.


---

## Conclusions

The differential expression analysis of ESR1 wild-type, Y537S, and D538G mutant T47D breast cancer cells reveals that both mutations cause extensive transcriptional reprogramming, affecting over 60% of the expressed transcriptome. The Y537S mutation produces a slightly broader transcriptional impact (13,914 DE genes) compared to D538G (13,149 DE genes), consistent with previous reports suggesting that Y537S may exhibit a more distinct and aggressive molecular phenotype. Both mutations appear to constitutively influence classical estrogen target genes involved in cell proliferation, survival, and breast cancer progression, even in the absence of estrogen.

The treatment time comparison identifies a smaller but relevant set of 5,614 time-responsive genes, indicating that mutation status is likely the primary driver of transcriptomic variation, while temporal estrogen signaling may act as a secondary modulatory factor. The heatmap analysis shows that ESR1 mutant cells cluster distinctly from wild-type cells, suggesting a fundamental shift in gene regulatory programs associated with these clinically relevant mutations.

Together, these findings support the possibility that ESR1 ligand-binding domain mutations Y537S and D538G contribute to endocrine therapy resistance in metastatic ER-positive breast cancer through sustained, ligand-independent transcriptional activity. This transcriptional rewiring may have important implications for the development of therapeutic strategies aimed at targeting mutant-specific estrogen receptor signaling.
