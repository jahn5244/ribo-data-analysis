
# Ribo-Data Analysis with NOISeq

This repository contains an analysis pipeline for evaluating differential gene expression using the Ribo-ITP approach and the NOISeq R package.

## Project Context

**Objective**: To investigate how gene translation changes across four timepoints following Long-Term Potentiation (LTP), a neural process tied to memory and learning.

**Technology Used**:
- Ribo-ITP: Novel ribosome profiling using isotachophoresis
- NOISeq: RNA-seq quality control and DEG detection
- R / Bioconductor
- Ensembl annotations

## Dataset

- 4 timepoints: Control, 10 min, 30 min, 1 hr
- 3 biological replicates: A, B, C (only A retained after QC)
- Data source: [GEO GSE185732](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185732)

---

## 🔬 Results Summary

### ✅ Replicate Filtering
Only **replicate A** was retained based on read length distribution. Ideal ribosome footprints are ~28 nt. B and C replicates showed irregular profiles.

![Replicate Filtering](/images/page6_img1.jpeg)

### 📉 Jaccard Index Between Time Points

| Comparison        | Jaccard Index |
|-------------------|---------------|
| Control vs 10 min | 0.45          |
| 10 min vs 30 min  | 0.13          |
| 30 min vs 1 hr    | 0.30          |
| Control vs 1 hr   | 0.36          |

🧠 **Interpretation**: Large shift between 10 → 30 min after LTP indicates an important regulatory transition.

### 📊 RNA-Seq QC with NOISeq

#### Biotype Detection
- Most transcripts were protein-coding
- Low contamination risk

![Biotype Detection](/images/page11_img1.jpeg)
![More QC](/images/page11_img2.jpeg)

#### Sequencing Bias Checks
- **Length bias**: Slight correlation (high R²)
- **GC content bias**: Minimal
- **RNA composition bias**: Detected (FAILED diagnostic test)

![Bias Detection - Length](/images/page12_img1.jpeg)
![Bias Detection - GC Content](/images/page12_img2.jpeg)

### 📈 Differential Expression Results

| Comparison        | DEGs |
|-------------------|------|
| Control vs 10 min | 854  |
| 10 min vs 30 min  | 331  |
| 30 min vs 1 hr    | 786  |
| Control vs 1 hr   | 702  |

#### 🔼 Upregulated / 🔽 Downregulated Summary

| Comparison        | 🔽 Down | 🔼 Up | Ratio (Up/Down) |
|-------------------|--------|-------|------------------|
| Control vs 10 min | 684    | 170   | 0.25             |
| 10 min vs 30 min  | 85     | 246   | 2.89             |
| 30 min vs 1 hr    | 672    | 114   | 0.17             |
| Control vs 1 hr   | 679    | 23    | 0.03             |

---

## 📚 References

- [NOISeq Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666377/)
- [Ribo-ITP Nature Paper](https://www.nature.com/articles/s41586-023-06228-9)
- [Dataset GEO: GSE185732](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185732)

---
