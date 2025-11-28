# Gene Expression Analysis of Retinal Tissues in Glaucoma
**A Bioinformatics (RNA-seq analysis) & Machine Learning Pipeline Using Public RNA-seq Data**

---

## Project Overview

This repository contains a portfolio-grade, end-to-end RNA-seq + machine learning pipeline **(QC → quantification → DE → pathway analysis → ML)** built around a public mouse retina dataset (GEO: **GSE212186**) that models:
- **Physiological aging** of the retina, and  
- **Glaucoma-like stress** via unilateral intraocular pressure (IOP) elevation.

**Design**: This project is designed to be implemented in two runs of DESeq2 and ML models to generate a clear biological story that explains our dataset clearly.
- **First Run**: Glaucoma (IOP-treated) vs. Untreated control retinas.
  
  This run focuses specifically on practicing and demonstrating a complete RNA-seq analysis pipeline:
  
  >- **Raw data → QC → Salmon → tximport → DESeq2 → pathway/TF scores → ML (logistic, RF) → nested CV.**
  >- **Differential expression (DE) and classification of IntraOcular Pressure (IOP)–elevated retinas (glaucoma) vs contralateral control retinas.**
  >- **Derive **biologically meaningful features** (pathway / TF / signalling activities) from the gene expression data.**
  >- **Train **machine learning classifiers** (regularized logistic regression and random forest) to distinguish **IOP vs control** using **nested cross-validation**.**

- **Final Run**: Focuses specifically on clearly separating samples by group to find Glaucoma-like stress effect while aging (young vs old) and IOP stress elevation effect at each age separately.

**Final Goal**:
- Adjust our metadata into different groups other than only IOP (glaucoma) vs control to tell a clear biological story:
  
  >- **Physiological aging effect (old vs young) in retina.**
  >- **Glaucoma-like stress effect (IOP vs control) separately in young and old.**

All steps are implemented with the mindset of an educational, **portfolio-garde** project: clear directory structure and scripts; documented analysis & decisions.

---

## Dataset Description

- **Accession**: [GSE212186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212186)
- **Species**: *Mus musculus* (mouse)
- **Tissue**: Retina
- **Samples**: 20 RNA-seq runs from ~10 BioSamples (some are technical replicates)

Mouse retina bulk RNA-seq with unilateral IOP elevation (constant 30 mmHg for 1 hr) and matched untreated controls: Young/old cohorts noted in the paper description (we focused the first run on IOP vs. control).

---

## Project Pipeline

| Step                     | Tools Used                              |
|--------------------------|------------------------------------------|
| Metadata Curation        | Manual, `SraRunTable.csv`                |
| FASTQ Download           | `fasterq-dump`, `pigz`                   |
| Quality Control          | `FastQC`, `MultiQC`                      |
| Read Alignment           | Salmon (quasi-mapping/ pseudoalignment)                    |
| Quantification           | tximport                          |
| Differential Expression  | `DESeq2`                             |
| Functional Enrichment    | `clusterProfiler`         |
| Machine Learning         | glmnet (Regularized logistic with α–λ tuning, standardized), ranger (Random Forest with permutation importance)    |

---

## Download & QC

  -	Fetch reads (e.g., SRA Toolkit → FASTQ) into data/raw_fastq/.

---

## Metadata Curation

**First Run**:
  - samples.txt created with a sample column (sample IDs) and a condition column (IOP, control) from SraRunTable.csv.
  - We mapped any examination or treatment of “iop/treated” → IOP, “control/untreated/ctrl” → control.

---

## QC → Trimming → QC

  - Raw FASTQ files (20 RNA-seq runs, paired-end) were assessed using FastQC and MultiQC.
 	  - **FastQC:** per-sample quality metrics
    - **MultiQC:** aggregates across all samples
      
  - **Output:** HTML reports, per-base quality, GC content, adapter content
  - Overall high per-base quality (Phred > 30 across most of each read), with a mild drop at the 3′ ends typical for Illumina sequencing.
  - Noticeable adapter contamination and several overrepresented sequences consistent with Illumina adapters and poly(A/T) tails.
  - Systematic per-base sequence content bias in the first ~10 bp, likely due to random hexamer priming and residual adapter/primer sequence.
  - Elevated duplication levels, as expected in bulk RNA-seq from a complex tissue. Duplication was interpreted cautiously as a mixture of technical and biological duplication rather than an automatic QC failure.

**Trimming**: To address these issues, I applied read trimming using fastp with the following strategy:
  - Hard-trimmed the first 10 bases from both read 1 and read 2 to remove biased initial sequence content.
  - Enabled automatic adapter detection and removal for paired-end data.
  - Performed 3′-end quality trimming, removing low-quality bases (mean quality < 20).
  - Discarded reads shorter than 50 bp after trimming.

**After trimming, QC metrics improved**: Trimming was successful. Adapter contamination, low-quality tails, and abnormal GC content were effectively removed. This will improve the accuracy of downstream steps (alignment, quantification, DE analysis).

All QC reports (FastQC, MultiQC, and fastp HTML reports) are stored in results/qc/.

---

## Alignment + Quantification (Salmon)

  - We switched from STAR/featureCounts → Salmon quasi-mapping for pseudoalignment & quantification
  - **Goal**: Get robust gene-level abundances fast so we can move on to DE and pathway analysis. Quantify transcript abundance using pseudoalignment against a reference transcriptome.
.
  - **Why Salmon fits?**:
    - Pseudo mapping: no full alignment since no heavy BAMs
    -	Handles many biases (GC, positional, sequence): --validateMappings, --seqBias, --gcBias. 
    -	No need to deduplicate—counts are model-based, not read-by-read. 
    -	Produces transcript-level TPM and estimated counts, necessary for downstream differential expression.
    -	Plays perfectly with tximport/DESeq2 (industry standard).

  -	Index reference transcriptome (Gencode/Ensembl GRCm39 / mm39).
  -	Quantify each sample with Salmon (libtype auto-detect; validate mapping rates).
  - Salmon produced transcript-level quantification files (quant.sf), each containing: Name, Length, EffectiveLength, TPM, and NumReads.

All Outputs under results/salmon/.

---

## Import & DE (first run: IOP vs control)

- Each gene can have multiple transcript isoforms. However, most biological questions (like glaucoma vs control differences) are at the gene level. Hence, we aggregate **transcript abundances per gene** using:
  -	**tximport** (R package) → Reads in all quant.sf files.
  -	**Tx2Gene** mapping file → Maps transcript IDs → gene IDs.

- Now we have gene-level counts, ready for **statistical modeling**. Raw counts differ due to:
  - Sequencing depth (some samples have more reads)
  - RNA composition (some samples overexpress a few genes)
  - Gene length differences.

- **Normalization**: makes counts comparable across samples. Save normalized counts/variance-stabilized data for downstream.
- **DESeq2**:
  -  models read counts using a Negative Binomial distribution, which handles overdispersion (biological variability).
  -  estimates fold changes and p-values to find differentially expressed (DE) genes between groups.
  - **Outputs**:
    - Plots: PCA, MA, Volcano, Heatmaps
    - DEG table:  a list of Differentially Expressed Genes (DEGs) with log₂ fold change (log₂FC), p-values, and adjusted p-values (FDR) saved in a CSV file.

 All Outputs under results/deseq2/.
 
---

## Comprehensive Functional Enrichment Analysis: Gene Ontology and Pathways 

- After DESeq2, we have a list of genes with log₂ fold change and adjusted p-values (FDR). We'll use these for:
  -	**Gene Ontology (GO) enrichment**: Over-Representation Analysis (ORA) of GO terms, specifically for the three main categories:
    -	Biological Process (BP),
    -	Molecular Function (MF),
    -	Cellular Component (CC)
      
  -	**Pathway Analysis**: Reactome, KEGG.
  -	Biological interpretation (e.g., oxidative stress pathways in glaucoma).

This is where **data science meets biology**.

---

## Pathway / TF / signalling activity scoring

- From VST-normalized expression, we compute activity tables:

  - **DoRothEA:** TF activity scores per sample.  
  - **GSVA:** Hallmark gene set activity scores per sample.  
  - **PROGENy:** Signalling pathway activity scores per sample.
    
Each resulting file is a matrix of features × samples and is used as input to the ML step.

---

## Machine Learning  

- **First Run**: IOP (Glaucoma) vs control classification
  - Use the activity matrices as input features to classify **IOP vs control** retinas.
  - Models: Elastic-net logistic regression with α–λ tuning (glmnet), Random Forest with permutation importance (ranger).
  - Repeated CV (5×5): Pooled CV curves + 95% CI for model evaluation and feature importance.
  -	**nested cross-validation**: Unbiased generalization
    -	Hyperparameters tuning (5x5folds inner CV).
    - Unbiased estimate of predictive performance (5-fold outer CV).
    - Each block includes (i) **nested outer ROC**, (ii) **fold AUCs**, (iii) **mean±sd AUC**, (iv) **best inner-CV tunes**, (v) **outer predictions**.
      
  -	Feature stability (RF) via repeated subsampling.

---

## Biological Interpreation (Top Features)

### DoRothEA — Stable TFs
- **E2F4, IRF1, HIF1A, TP53, MITF, USF1** (frequently top-ranked across resamples).  
 
### GSVA (Hallmarks) — Stable Hallmark Pathways
- **Angiogenesis, EMT, Xenobiotic metabolism, Apoptosis, Cholesterol homeostasis, ROS, PI3K/AKT/mTOR, UV response (UP)**.  

### PROGENy — Stable Pathways
- **JAK/STAT, EGFR, Estrogen, TGF-β, WNT, MAPK, NFκB, Hypoxia**. 

---

## Key Results

```markdown

Feature space	Model	Repeated-CV AUC	Nested-CV AUC	Artifacts
DoRothEA (TF activity)	glmnet	0.809	0.830	ROC (repeated), ROC (nested), 95% CI ROC
RF	0.823	0.825	ROC (repeated), ROC (nested), 95% CI ROC
GSVA (MSigDB Hallmarks)	glmnet	0.833	0.780	ROC (repeated), ROC (nested), 95% CI ROC
RF	0.918	0.910	ROC (repeated), ROC (nested), 95% CI ROC
PROGENy (pathway activities)	glmnet	0.501	0.360	ROC (repeated), ROC (nested), 95% CI ROC
RF	0.841	0.915	ROC (repeated), ROC (nested), 95% CI ROC
GSVA + PROGENy (meta-RF)	RF ensemble	0.921 (95% CI 0.900–0.943)	—	ROC, summary
```



```
Feature space	Model	Repeated-CV AUC	Nested-CV AUC	Artifacts
DoRothEA (TF activity)	glmnet	0.809	0.830	ROC (repeated), ROC (nested), 95% CI ROC
	RF	0.823	0.825	ROC (repeated), ROC (nested), 95% CI ROC
GSVA (MSigDB Hallmarks)	glmnet	0.833	0.780	ROC (repeated), ROC (nested), 95% CI ROC
	RF	0.918	0.910	ROC (repeated), ROC (nested), 95% CI ROC
PROGENy (pathway activities)	glmnet	0.501	0.360	ROC (repeated), ROC (nested), 95% CI ROC
	RF	0.841	0.915	ROC (repeated), ROC (nested), 95% CI ROC
GSVA + PROGENy (meta-RF)	RF ensemble	0.921 (95% CI 0.900–0.943)	—	ROC, summary
Takeaways. Nested-CV confirms strong generalization for GSVA-RF (~0.91 AUC) and PROGENy-RF (~0.92 AUC). Linear models underperform on PROGENy → signal is non-linear. A simple RF ensemble of GSVA+PROGENy reaches AUC ≈ 0.921.

```

---

## Repository structure
```text
glaucoma_rnaseq_pipeline/
├── data/
│   ├── raw_fastq/            # FASTQ or SRA-derived reads
│   ├── metadata/             # samples.txt, SRP394552_metadata.csv, metadata_deseq2.csv metadata files
│   └── trimmed_fastq/        # FASTQ reads after trimming
│   └── reference/            # tx2gene.csv maps for tximport, downloaded reference genome, built index
├── results/
│   ├── qc/                   # FastQC / MultiQC reports
│   ├── salmon/               # Salmon quant outputs per sample, post-salmon QC results
│   ├── deseq2/               # DE analysis: DEGs table and PCA plot
│   ├── pathways/             # GO enrichment analysis results( ORA GO: BP, MF, CC); Pathway analysis results (Reactome, KEGG); dot plots, emaps; DoRothEA TF activity, GSVA Hallmark, PROGENy
│   └── delivery/             # Export enrichment analysis matrices (counts, normalized, VST, TPM); DE with covariate + LFC shrinkage results
│   └── ml/                   # ML model (glmnet, RF) results
│   └── ml_eval/              # Repeated 5x5 CV ML results
│   └── nested_cv/            # nested CV results
├── scripts/                  # python, bash, and R scripts
└── README.md
```

---

## Clone the Repo
```bash
git clone https://github.com/yourusername/glaucoma-rnaseq-pipeline.git
cd glaucoma-rnaseq-pipeline
```
---

## Acknowledgement

Used AI tools (ChatGPT & Gemini) to plan & design project; improve & test code; prepare & refine readme.

---

## 🧠 Author
Md Ishtiaque Hossain \
MSc Candidate, Computer and Information Sciences \
University of Delaware \
[LinkedIn](linkedin.com/in/ishtiaque-h) | [GitHub](https://github.com/Ishtiaque-h)
