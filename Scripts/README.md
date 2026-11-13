# Analysis of Inflammatory and Immune-Related Gene Expression Changes Following Endurance Exercise Using Bioinformatics

## 1. Introduction

Physical exercise produces several physiological and molecular responses in the human body. Endurance exercise, in particular, can influence immune-system activity, inflammatory responses and cellular metabolism. These responses involve changes in the expression of genes associated with different biological processes.

Gene expression analysis provides an opportunity to investigate how exercise influences molecular activities by examining changes in RNA expression levels. Bioinformatics techniques can be used to identify differentially expressed genes, investigate their biological functions and identify associated molecular pathways.

This project investigates inflammatory and immune-related gene expression changes following endurance exercise using publicly available gene expression data from the Gene Expression Omnibus (GEO) database. The analysis focuses on circulating neutrophil samples collected at different time points following exercise.

The project applies several bioinformatics techniques, including exploratory data analysis, Principal Component Analysis (PCA), sample correlation analysis, differential gene expression analysis, gene annotation, Gene Ontology (GO) enrichment and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway analysis.

## 2. Research Objectives

The main objective of this project is to investigate inflammatory and immune-related gene expression changes following endurance exercise using bioinformatics approaches.

The specific objectives are:

- To obtain and prepare publicly available gene expression data related to endurance exercise.
- To investigate gene expression patterns using exploratory data analysis.
- To identify differentially expressed genes at different post-exercise time points.
- To annotate differentially expressed probes with corresponding gene information.
- To identify enriched biological processes using GO enrichment analysis.
- To identify associated biological pathways using KEGG pathway analysis.
- To investigate selected immune-related genes and compare their expression patterns across different post-exercise time points.

## 3. Dataset Description

The dataset used in this project was obtained from the Gene Expression Omnibus (GEO) database.

**Dataset accession:** GSE43856

**Platform:** GPL10558 – Illumina HumanHT-12 V4.0 Expression BeadChip

**Organism:** Homo sapiens

The original study investigated transcriptome responses to endurance exercise in healthy endurance-trained men. The exercise protocol involved one hour of cycling followed by one hour of running.

The study included eight healthy endurance-trained male participants. Blood polymorphonuclear neutrophil (PMN) samples and skeletal muscle samples were collected at different time points.

For this project, the analysis focused on the circulating blood PMN samples collected at four time points:

- Before exercise (baseline)
- 3 hours after exercise
- 48 hours after exercise
- 96 hours after exercise

This resulted in 32 blood PMN samples for the selected analysis.

The dataset provides an opportunity to investigate changes in gene expression following endurance exercise and examine how selected gene expression patterns vary over time.

## 4. Methodology

The analysis was performed using the R programming language and RStudio. Several R packages were used for data processing, statistical analysis, gene annotation and biological interpretation.

### 4.1 Data Preparation

The gene expression dataset and corresponding platform annotation information were obtained from GEO.

The expression matrix was prepared by selecting the 32 blood PMN samples relevant to the research objectives. Sample information was examined to identify the corresponding collection time points.

The GEO expression matrix had already undergone background correction and quantile normalization as part of the original data processing. Therefore, additional normalization was not applied unnecessarily.

### 4.2 Exploratory Data Analysis

Exploratory data analysis was performed to understand the distribution and characteristics of the gene expression values.

Expression distribution plots were generated to examine the overall patterns in the dataset. These visualisations helped assess the consistency of expression values across the selected samples.

### 4.3 Principal Component Analysis

Principal Component Analysis (PCA) was performed to investigate similarities and differences among the samples.

PCA reduces the dimensionality of high-dimensional gene expression data while preserving important variation. The resulting plot provides a visual representation of sample relationships.

This analysis helped explore whether samples collected at different time points showed differences in their overall gene expression patterns.

### 4.4 Sample Correlation Analysis

Sample correlation analysis was performed to examine similarities between the selected samples.

A correlation heatmap was generated to visualise the relationships between samples. Highly correlated samples generally have more similar gene expression patterns.

This analysis provided additional information about sample consistency and overall expression relationships.

### 4.5 Differential Gene Expression Analysis

Differential gene expression analysis was performed to identify genes with changes in expression following endurance exercise.

Each post-exercise time point was compared with the baseline condition:

- 3 hours versus baseline
- 48 hours versus baseline
- 96 hours versus baseline

The analysis used the limma statistical framework for microarray gene expression analysis.

Log-fold change values were used to describe the direction and magnitude of expression differences. Multiple-testing correction was applied using the Benjamini-Hochberg method.

Genes with adjusted p-values below 0.05 were considered statistically significant.

Volcano plots were generated to visualise the relationship between statistical significance and the magnitude of gene expression changes.

### 4.6 Gene Annotation

The differentially expressed probes were annotated using the platform information.

Probe identifiers were mapped to corresponding gene symbols, Entrez Gene IDs and gene descriptions.

This step helped convert probe-level results into biologically interpretable gene information.

### 4.7 Gene Ontology Enrichment Analysis

Gene Ontology (GO) enrichment analysis was performed using the `clusterProfiler` package and the `org.Hs.eg.db` annotation database.

The analysis focused on the Biological Process (BP) category to identify biological processes overrepresented among differentially expressed genes.

Separate enrichment analyses were also performed for upregulated and downregulated genes.

The Benjamini-Hochberg method was used to adjust p-values. Enriched terms were examined to understand the biological processes associated with the observed expression changes.

### 4.8 KEGG Pathway Enrichment Analysis

KEGG pathway enrichment analysis was performed to identify biological pathways associated with the differentially expressed genes.

The analysis used the `enrichKEGG` function from the `clusterProfiler` package with the human organism identifier.

The identified pathways were examined to understand the molecular processes represented among the differentially expressed genes.

### 4.9 Immune-Related Gene Analysis

Ten immune-related genes were selected for further investigation based on the identified Th1 and Th2 cell differentiation pathway.

The selected genes were:

- RUNX3
- CD247
- CD3D
- STAT4
- NFATC1
- PPP3CC
- IL2RB
- MAPK14
- HLA-DRA
- HLA-DMB

Their expression changes were investigated at 3 hours after exercise. The selected genes were also compared across the 3-hour, 48-hour and 96-hour time points.

For the time-point comparison, median log-fold change values across probes mapped to each gene were used as descriptive summaries.

## 5. Results and Findings

### 5.1 Exploratory Data Analysis

Expression distribution plots were generated to examine the gene expression values across the selected samples.

PCA and sample correlation heatmaps were also generated to explore similarities and differences among samples.

These analyses provided an overview of the expression data before performing differential expression and enrichment analyses.

### 5.2 Differential Gene Expression

Differential expression analysis identified 565 significant probe-level results at 3 hours after exercise using an adjusted p-value threshold of 0.05.

Among these results:

- 44 probe rows showed positive log-fold changes.
- 521 probe rows showed negative log-fold changes.

These values represent significant probe-level results rather than unique genes. Some genes may be represented by multiple probes.

The results were annotated to identify corresponding gene symbols and gene identifiers.

Separate differential expression analyses were also performed for the 48-hour and 96-hour comparisons.

### 5.3 GO Enrichment Results

The GO Biological Process enrichment analysis identified 50 significantly enriched terms among the differentially expressed genes at 3 hours.

Some of the most significant terms included:

- T-cell activation
- Alpha-beta T-cell activation
- tRNA metabolic process
- tRNA processing
- Mitochondrial gene expression
- RNA modification
- Ribosome biogenesis
- Adaptive immune response

The separate analysis of upregulated and downregulated genes identified 45 enriched terms among downregulated genes. No significantly enriched GO terms were identified among upregulated genes.

The results indicate that several immune-related and cellular processes were represented among the differentially expressed genes.

### 5.4 KEGG Pathway Results

KEGG enrichment analysis identified eight significantly enriched pathways at 3 hours after exercise.

The identified pathways were:

1. Ribosome biogenesis in eukaryotes
2. Base excision repair
3. Prion disease
4. Amyotrophic lateral sclerosis
5. Th1 and Th2 cell differentiation
6. Aminoacyl-tRNA biosynthesis
7. Spliceosome
8. Carbon metabolism

These pathways represent several biological areas, including RNA processing, protein synthesis, DNA repair, metabolism and immune-related functions.

The identification of disease-named pathways does not indicate that participants had those diseases. Such pathways can share genes involved in general cellular processes.

### 5.5 Immune-Related Gene Expression Results

The expression analysis of the ten selected immune-related genes showed statistically significant differential expression at 3 hours after exercise.

Eight selected genes showed negative log-fold changes:

- RUNX3
- CD247
- CD3D
- STAT4
- PPP3CC
- IL2RB
- HLA-DRA
- HLA-DMB

Two genes showed positive log-fold changes:

- NFATC1
- MAPK14

These findings demonstrate that the selected genes exhibited different expression patterns following endurance exercise.

### 5.6 Comparison Across Post-Exercise Time Points

The expression patterns of the selected immune-related genes were compared across 3, 48 and 96 hours after exercise.

At 3 hours, all ten selected genes showed statistically significant differential expression based on their selected representative probe results.

However, none of the selected genes showed statistically significant differential expression at 48 or 96 hours after multiple-testing correction.

Several genes also showed log-fold changes closer to zero at the later time points.

This suggests that some early expression changes may have moved closer to baseline over time. However, the analysis did not directly test differences between post-exercise time points. Therefore, statistical recovery or return to baseline cannot be established from these comparisons alone.

## 6. Biological Interpretation and Discussion

The analysis identified changes in the expression of several genes in circulating neutrophils following endurance exercise.

The GO enrichment results highlighted immune-related annotations, including T-cell activation and adaptive immune response. Other enriched terms involved tRNA metabolism, RNA processing, mitochondrial gene expression and ribosome biogenesis.

These findings suggest that the observed gene expression changes were associated with several immune-related and cellular processes.

The KEGG results provided additional information about the pathways represented among the differentially expressed genes. In particular, the enrichment of the Th1 and Th2 cell differentiation pathway highlighted the presence of immune-related genes in the results.

The selected immune-related genes showed significant expression changes at 3 hours after exercise. However, it is important to consider the biological context of the dataset. The analysed samples were circulating neutrophils rather than isolated T cells. Therefore, enrichment of a T-cell-associated pathway does not directly demonstrate changes in T-cell differentiation or activity.

Furthermore, GO and KEGG enrichment analyses identify overrepresented annotations and pathways. They do not directly measure whether biological processes are activated or suppressed.

The time-point comparison showed that none of the selected immune-related genes remained statistically significant at 48 or 96 hours compared with baseline. Although some expression changes appeared closer to zero at later time points, further statistical investigation would be necessary to establish temporal differences.

Overall, the results provide evidence of gene expression responses associated with endurance exercise, particularly at the early post-exercise time point.

## 7. Limitations

Several limitations should be considered when interpreting the findings.

1. **Small sample size:** The study included eight participants, which limits the statistical power and generalisability of the findings.

2. **Specific sample type:** The analysis focused only on circulating neutrophil samples. Therefore, the findings cannot represent the responses of all immune cell types.

3. **Probe-level analysis:** Multiple probes may map to the same gene. Consequently, probe-level results do not always represent unique genes.

4. **Gene expression measurements:** Changes in RNA expression do not necessarily indicate corresponding changes in protein levels or biological activity.

5. **Pathway interpretation:** Enrichment results indicate overrepresentation of genes associated with biological annotations rather than direct evidence of pathway activation or suppression.

6. **Time-point comparisons:** The comparisons were performed against baseline separately. Direct statistical tests between post-exercise time points would be required to establish significant temporal differences.

7. **Selected genes:** The detailed immune-related analysis focused on ten selected genes and does not represent every inflammatory or immune-related gene in the dataset.

## 8. Conclusion

This project investigated inflammatory and immune-related gene expression changes following endurance exercise using the publicly available GSE43856 dataset.

The analysis included data preparation, exploratory data analysis, PCA, sample correlation, differential gene expression analysis, gene annotation, GO enrichment and KEGG pathway enrichment.

The results identified 565 significant probe-level results at 3 hours after exercise. GO enrichment analysis identified 50 enriched biological processes, while KEGG analysis identified eight enriched pathways.

Further investigation of ten selected immune-related genes showed significant expression changes at 3 hours. However, none of these selected genes showed statistically significant differential expression at 48 or 96 hours compared with baseline.

The findings demonstrate that endurance exercise was associated with changes in the expression of several genes in circulating neutrophils. They also highlight the usefulness of bioinformatics techniques for investigating molecular responses to physical exercise.

Future studies involving larger participant groups, additional immune cell types and functional measurements could provide a more comprehensive understanding of exercise-related immune responses.

## 9. Tools and Technologies

The following tools and technologies were used:

- **R:** Statistical computing and data analysis.
- **RStudio:** Development environment for R scripts.
- **GEO:** Source of the gene expression dataset.
- **limma:** Differential gene expression analysis.
- **clusterProfiler:** GO and KEGG enrichment analysis.
- **org.Hs.eg.db:** Human gene annotation.
- **ggplot2:** Data visualisation.
- **pheatmap:** Sample correlation heatmap.
- **Git and GitHub:** Version control and project management.

## 10. Repository Structure

```text
2021-ASP-38/
│
├── Scripts/
│   ├── 1_data_preparation.R
│   ├── 02_PCA_analysis.R
│   ├── 03_sample_correlation.R
│   ├── 04_metadata_check.R
│   ├── 05_differential_expression.R
│   ├── 06_probe_gene_annotation.R
│   ├── 07_GO_enrichment.R
│   ├── 08_KEGG_enrichment.R
│   └── 09_immune_gene_time_comparison.R
│
├── data/
│   └── raw/
│
├── figures/
│
├── results/
│
└── README.md
```

## 11. References

1. Neubauer, O., König, D., and Wagner, K. H. (2013). Recovery after an Ironman triathlon: sustained inflammatory responses and haemolysis. *European Journal of Applied Physiology*. The original study associated with GEO dataset GSE43856 should be checked against the GEO record before finalising the complete bibliographic details.

2. Gene Expression Omnibus (GEO). Dataset GSE43856. National Center for Biotechnology Information (NCBI).

3. Ritchie, M. E., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47.

4. Yu, G., Wang, L. G., Han, Y., and He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS: A Journal of Integrative Biology*, 16(5), 284–287.