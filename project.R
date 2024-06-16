install.packages("BiocManager")
BiocManager::install(c("ribor", "NOISeq", "biomaRt"))
install.packages("tidyverse")
devtools::install_version("dbplyr", version = "2.3.4")
library(ribor)
library(NOISeq)
library(biomaRt)
library(tidyverse)

# Load ribo file
file.path <- "/Users/johnahnline/Downloads/all.ribo" # Change based on file location
original.ribo <- Ribo(file.path, rename = rename_default)

# Plot length distribution
plot_length_distribution(x = original.ribo,
                         region = "CDS",
                         range.lower = 15, 
                         range.upper = 40,
                         fraction = TRUE)

# Plot region counts for quality control
plot_region_counts(x   = original.ribo,
                   range.lower = 25,
                   range.upper = 31)


# Prepare RNA-seq data 
experiments <- c("WT_control_A", "WT_10min_A", "WT_30min_A", "WT_1hr_A")
rnaseq <- get_rnaseq(ribo.object = original.ribo,
                     tidy = FALSE,
                     alias = TRUE,
                     experiment = experiments)

rnaseq <- rnaseq[, c("transcript", "experiment", "CDS")]
mycounts <- as.data.frame(rnaseq)
mycounts <- mycounts %>% # Rearrange mycounts
  pivot_wider(names_from = experiment, values_from = CDS, values_fill = 0) 

transcripts_order <- mycounts$transcript # To rearrange other DataFrames

# Load transcript data
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") 
transcript_data <- getBM(attributes = c("external_transcript_name", 
                                        "gene_biotype", "percentage_gene_gc_content"),
                         mart = ensembl)

# Find biotypes
biotypes <- transcript_data %>%
  select(external_transcript_name, gene_biotype) %>%
  distinct()
biotypes <- biotypes[!duplicated(biotypes$external_transcript_name), ]

transcript_names <- rownames(mycounts)
biotypes_filtered <- biotypes[biotypes$external_transcript_name %in% transcript_names, ]
order_indices <- match(rownames(mycounts), biotypes_filtered$external_transcript_name)
biotypes <- biotypes_filtered[order_indices, ]

# Find GC content
gc_content <- transcript_data %>%
  select(external_transcript_name, percentage_gene_gc_content) %>%
  distinct()
gc_content <- gc_content[!duplicated(gc_content$external_transcript_name), ]

gc_filtered <- gc_content[gc_content$external_transcript_name %in% transcript_names, ]
order_indices <- match(rownames(mycounts), gc_filtered$external_transcript_name)
gc_content <- gc_filtered[order_indices, ]

# Find CDS lengths of transcripts for NOISeq
region_coord <- get_internal_region_coordinates(ribo.object = original.ribo, alias = TRUE)
mylength <- region_coord %>%
  mutate(CDS_length = CDS_stop - CDS_start + 1) %>%
  select(transcript, CDS_length)
mylength<- mylength %>%
  arrange(factor(transcript, levels = transcripts_order))

# Organize factors for NOISeq
# From NOISeq documentation:
# Variables indicating the experimental group for each sample
# Must have as many rows as samples (columns in data object) = 4
myfactors = data.frame(Time = c("control", "10min", "30min", "1hr"),
                     Replicate = c("A", "A", "A", "A"))

# NOISeq object
mycounts <- mycounts %>% # Set transcript as index
  tibble::column_to_rownames(var = "transcript")
mydata <- readData(data = mycounts, factors = myfactors, length = mylength)
                  #biotype = biotypes, gc = gc_content) # ENSEMBL DOWN

#### Quality control

## Biotypes
# Biodetection plot
biodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
explo.plot(biodetection, samples = 4, plottype = "persample") # Change samples
explo.plot(biodetection, samples = c(1, 2), toplot = "protein_coding", 
           plottype = "comparison")
# Count distribution
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = c(1, 2, 3, 4), plottype = "boxplot")

##Sequencing Depth & Expression Quantification
# Saturation plot
saturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(saturation, toplot = 1, samples = 1:4, yleftlim = NULL, yrightlim = NULL)
# Count distribution per sample
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
# Sensitivity plot
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")






### Filtration

mycountsfilt = filtered.data(mycounts, factor = myfactors$Time, norm = FALSE,
                             depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")

transcriptsfilt<- rownames(mycountsfilt)
mylengthfilt <- mylength %>%
  filter(transcript %in% transcriptsfilt)
mydatafilt <- readData(data = mycountsfilt, factors = myfactors, length = mylengthfilt)




### Sequencing bias detection
# Length bias
mylengthbias = dat(mydatafilt, factor = "Replicate", type = "lengthbias")
explo.plot(mylengthbias, toplot = "global")

myRPKM = tmm(assayData(mydatafilt)$exprs, long = mylengthfilt, lc = 0)
mydatacorr <- readData(data = myRPKM, factors = myfactors, 
                       length = mylengthfilt)
mylengthbiascorr = dat(mydatacorr, factor = "Replicate", type = "lengthbias")
explo.plot(mylengthbiascorr, toplot = "global")  


# GC content bias
myGCbias = dat(mydata, factor = "Replicate", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")
# RNA composition
mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd)

### Batch effect
# PCA
set.seed(123)
mycounts2 = mycounts
mycounts2 = mycounts2 + runif(nrow(mycounts2) * 4, 3, 5)
mydata2 = readData(mycounts2, factors = myfactors)
myPCA = dat(mydata2, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA, factor = "Time")

# Quality report
QCreport(mydata, samples = NULL, factor = "Time", norm = FALSE)



# NOISeq-Sim using raw data
# Recommendated parameters: replicates = "no", nss = 5, pnr = 0.2, v = 0.02

myconditions = c("control", "1hr") # Change
mynoiseq <- noiseq(mydata, factor = "Time", conditions = myconditions,
                   k = NULL, norm = "n", pnr = 0.2,
                   nss = 5, v = 0.02, lc = 0, replicates = "no") 
DE.plot(mynoiseq, q = 0.9, graphic = "expr", log.scale = TRUE)

myconditions = c("control", "1hr") # Change
mynoiseq <- noiseq(mydata, factor = "Time", conditions = myconditions,
                    k = NULL, norm = "n", pnr = 0.2,
                    nss = 5, v = 0.02, lc = 0, replicates = "no") 
DE.plot(mynoiseq, q = 0.9, graphic = "expr", log.scale = TRUE)


# Correcting data





### Normalization

myRPKM = rpkm(assayData(mydatafilt)$exprs, long = cds_lengthfilt, k = 0, lc = 1)
mydatacorr <- readData(data = myRPKM, factors = myfactors, 
                       length = cds_lengthfilt)



# NOISeq

myconditions = c("10min", "30min") # Change
mynoiseqcorr <- noiseq(mydatacorr, factor = "Time", conditions = myconditions,
                       k = NULL, norm = "n", pnr = 0.2,
                       nss = 5, v = 0.02, lc = 0, replicates = "no")
DE.plot(mynoiseqcorr, q = 0.8, graphic = "expr", log.scale = TRUE)

mynoiseq.deg = degenes(mynoiseqcorr, q = 0.8, M = NULL)



# Normalization
myRPKM1 = rpkm(assayData(mydata)$exprs, long = cds_length, k = 0, lc = 1)
mydatacorr1 <- readData(data = myRPKM1, factors = myfactors, 
                        length = cds_length)
# Length bias
mylengthbiascorr = dat(mydatacorr, factor = "Replicate", type = "lengthbias")
explo.plot(mylengthbiascorr, toplot = "global") 
# GC content bias
myGCbias = dat(mydata, factor = "Replicate", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")
# RNA composition
mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd)

myRPKM = uqua(assayData(mydata)$exprs, long = mylengthfilt, lc = 0)
mydatacorr <- readData(data = myRPKM, factors = myfactors, 
                       length = mylengthfilt)
mylengthbiascorr = dat(mydatacorr, factor = "Replicate", type = "lengthbias")
explo.plot(mylengthbiascorr, toplot = "global") 