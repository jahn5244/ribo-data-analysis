# Ribosome Profiling Data Analysis

## Project Overview

This R script is used for analyzing ribosome profiling and RNA-seq data to study gene expression and perform quality control checks. It performs various bioinformatics analyses including data preprocessing, quality control, and differential expression analysis using the `ribor`, `NOISeq`, and `biomaRt` packages.

The dataset can be accessed and downloaded using this link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185732 

## Setup Instructions

1. **Install Required Packages:**

    ```r
    install.packages("BiocManager")
    BiocManager::install(c("ribor", "NOISeq", "biomaRt"))
    install.packages("tidyverse")
    devtools::install_version("dbplyr", version = "2.3.4")
    ```

2. **Load Libraries:**

    ```r
    library(ribor)
    library(NOISeq)
    library(biomaRt)
    library(tidyverse)
    ```

## Data Loading and Preprocessing

- **Load Ribo Data:**

    ```r
    file.path <- "/path/to/all.ribo"  # Change based on file location for testing.
    original.ribo <- Ribo(file.path, rename = rename_default)
    ```

- **Plot Length Distribution:**

    ```r
    plot_length_distribution(x = original.ribo, region = "CDS", range.lower = 15, range.upper = 40, fraction = TRUE)
    ```

- **Quality Control:**

    ```r
    plot_region_counts(x = original.ribo, range.lower = 25, range.upper = 31)
    ```

## RNA-seq Data Preparation

- **Define Experiments and Load Data:**

    ```r
    experiments <- c("WT_control_A", "WT_10min_A", "WT_30min_A", "WT_1hr_A")
    rnaseq <- get_rnaseq(ribo.object = original.ribo, tidy = FALSE, alias = TRUE, experiment = experiments)
    rnaseq <- rnaseq[, c("transcript", "experiment", "CDS")]
    mycounts <- as.data.frame(rnaseq) %>%
      pivot_wider(names_from = experiment, values_from = CDS, values_fill = 0)
    ```

- **Load and Prepare Transcript Data:**

    ```r
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    transcript_data <- getBM(attributes = c("external_transcript_name", "gene_biotype", "percentage_gene_gc_content"), mart = ensembl)
    ```

## Data Analysis and Quality Control

- **Biotype and GC Content Filtering:**

    ```r
    biotypes <- transcript_data %>% select(external_transcript_name, gene_biotype) %>% distinct()
    gc_content <- transcript_data %>% select(external_transcript_name, percentage_gene_gc_content) %>% distinct()
    ```

- **Create NOISeq Object:**

    ```r
    myfactors <- data.frame(Time = c("control", "10min", "30min", "1hr"), Replicate = c("A", "A", "A", "A"))
    mycounts <- mycounts %>% tibble::column_to_rownames(var = "transcript")
    mydata <- readData(data = mycounts, factors = myfactors, length = mylength)
    ```

## Quality Control Plots

- **Biotype Detection:**

    ```r
    biodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
    explo.plot(biodetection, samples = 4, plottype = "persample")
    ```

- **Sequencing Depth and Expression Quantification:**

    ```r
    saturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
    explo.plot(saturation, toplot = 1, samples = 1:4)
    ```

## Filtration and Normalization

- **Filter Data:**

    ```r
    mycountsfilt <- filtered.data(mycounts, factor = myfactors$Time, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")
    mydatafilt <- readData(data = mycountsfilt, factors = myfactors, length = mylengthfilt)
    ```

- **Bias Detection:**

    ```r
    mylengthbias <- dat(mydatafilt, factor = "Replicate", type = "lengthbias")
    explo.plot(mylengthbias, toplot = "global")
    ```

- **GC Content Bias:**

    ```r
    myGCbias <- dat(mydata, factor = "Replicate", type = "GCbias")
    explo.plot(myGCbias, samples = NULL, toplot = "global")
    ```

## Differential Expression Analysis

- **Perform NOISeq Analysis:**

    ```r
    myconditions <- c("control", "1hr")
    mynoiseq <- noiseq(mydata, factor = "Time", conditions = myconditions, k = NULL, norm = "n", pnr = 0.2, nss = 5, v = 0.02, lc = 0, replicates = "no")
    DE.plot(mynoiseq, q = 0.9, graphic = "expr", log.scale = TRUE)
    ```

## Conclusion

- **Results and Further Steps:**

    - Continue with data normalization and additional bias correction.
    - Explore other analytical methods for improved accuracy.

## Requirements

- **R Packages:** `BiocManager`, `ribor`, `NOISeq`, `biomaRt`, `tidyverse`, `devtools`

## Usage

1. **Install Packages:**
    ```r
    install.packages("BiocManager")
    BiocManager::install(c("ribor", "NOISeq", "biomaRt"))
    install.packages("tidyverse")
    devtools::install_version("dbplyr", version = "2.3.4")
    ```

2. **Run Analysis Script:**

    ```r
    source("analysis_script.R")  # Update path to your script
    ```

## License

This project is licensed under the MIT License.
