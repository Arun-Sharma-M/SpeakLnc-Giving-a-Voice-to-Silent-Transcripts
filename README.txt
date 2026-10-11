### **README: 

**Author:** Arun Sharma M
**Version:** 1.0
**Last Updated:** October 2026

---

#### **Overview**

This document provides a comprehensive, step-by-step description of the computational and analytical workflow used for transcriptome-wide analysis of long non-coding RNAs (lncRNAs) and mRNAs in cholangiocarcinoma (CHOL). It details all major stages — from RNA-Seq data retrieval and preprocessing to differential expression, enrichment, survival, and network analyses.

---

#### **Contents**

The file is divided into the following main sections:

1. **Retrieval and Processing of RNA-Seq Data**

   * Acquisition of SRA and TCGA-CHOL datasets
   * Quality control using FastQC
   * Reference genome indexing and read alignment using HISAT2
   * BAM conversion, sorting, and read quantification via *samtools* and *bedtools*

2. **Global Transcriptome Analysis using DESeq2**

   * Differential expression analysis
   * Normalization, fold-change shrinkage, and statistical testing
   * Visualization: MA, Volcano, PCA, Mean-SD, and Heatmaps
   * Identification of top up/downregulated lncRNAs

3. **TPM Normalization and Fold Change Calculation**

   * Conversion of read counts to TPM values
   * Calculation of log2 fold changes and filtering (L2FC ≥ 2, p ≤ 0.055, TPM ≥ 1)
   * Generation of heatmap-ready datasets

4. **Screening of lncRNAs and Their mRNAs**

   * Comparative analysis between SRA and TCGA datasets
   * Venn diagram-based intersection of shared transcripts

5. **Survival Analysis**

   * Kaplan–Meier and log-rank tests for selected lncRNAs (e.g., *LINC02362*)
   * Integration of TCGA clinical metadata

6. **Functional and Pathway Enrichment**

   * GeneSetDB, KEGG, and Reactome analyses
   * Bubble plot generation for enriched biological processes

7. **RNA–RNA Interaction Analysis**

   * IntaRNA predictions for lncRNA–mRNA interactions (e.g., *lnc-ADH4-3:2* vs *ADH4*)
   * Visualization of minimal energy profiles

8. **Data Visualization and Validation**

   * Box-and-whisker plots of shortlisted transcripts with *p*-value annotations

---

#### **Software and Dependencies**

* **Operating System:** Ubuntu 24.04 LTS
* **Languages:** R (v4.4), Bash
* **Key R Packages:** `DESeq2`, `ggplot2`, `TCGAbiolinks`, `pheatmap`, `ComplexHeatmap`, `ggVennDiagram`, `survminer`, `impute`, `dplyr`, `tidyverse`
* **External Tools:** HISAT2, FastQC, samtools, bedtools, IntaRNA

---

#### **Usage**

This document serves as a reproducible methodological reference for bioinformatics workflows applied in the *SpeakLnc* project (*“Giving a Voice to Silent Transcripts”*).
Users can follow the stepwise R and shell scripts to reproduce each analytical stage, adapting parameters as needed for other RNA-Seq datasets.

---

