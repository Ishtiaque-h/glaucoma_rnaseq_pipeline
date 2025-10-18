
---

## 🧾 Dataset Description

- **Accession**: [GSE212186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212186)
- **Species**: *Mus musculus* (mouse)
- **Tissue**: Retina
- **Design**: Glaucoma (IOP-treated) vs. Untreated control
- **Samples**: 20 RNA-seq runs from ~10 BioSamples (some are technical replicates)

---

## ⚙️ Planned Pipeline

| Step                     | Tools Used                              |
|--------------------------|------------------------------------------|
| Metadata Curation        | Manual, `SraRunTable.csv`                |
| FASTQ Download           | `fasterq-dump`, `pigz`                   |
| Quality Control          | `FastQC`, `MultiQC`                      |
| Read Alignment           | `STAR` (or `HISAT2`)                     |
| Quantification           | `featureCounts`                          |
| Differential Expression  | `DESeq2` (R)                             |
| Functional Enrichment    | `clusterProfiler`, `g:Profiler`          |
| Machine Learning         | `scikit-learn`, `pandas`, `XGBoost`      |
| Reproducibility          | `venv`, optionally `Snakemake`           |

---

## 🚧 Current Progress

- [x] Project structure and plan finalized  
- [x] Metadata manually curated with group labels  
- [x] Download script for 20 FASTQ files created  
- [ ] Quality Control (FastQC + MultiQC)  
- [ ] Genome alignment with STAR  
- [ ] Count matrix generation  
- [ ] Differential expression analysis  
- [ ] Functional enrichment  
- [ ] ML classifier (RF, SVM, etc.)  
- [ ] Final visualization + dashboard (optional)

---

## 📦 Getting Started

### Clone the Repo

git clone https://github.com/yourusername/glaucoma-rnaseq-pipeline.git
cd glaucoma-rnaseq-pipeline


### Set Up Environment

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

### Download Data

bash scripts/download_fastq.sh data/metadata/GSE212186_metadata.csv

## 🧠 Author
Md Ishtiaque Hossain
MSc Candidate, Computer and Information Sciences
University of Delaware
LinkedIn | GitHub