# 2026 transcriptomics HMOs
Code for the analysis and plots of the transcriptomic data of D39 in 3 different medium (lactose, 2 HMOs).

# Data
The different conditions tested are defined in `data/conditions.tsv` and the samples in `data/samples.tsv`.

# RNA-seq data analysis
**Quality Control** : Reads QC was performed using [FastQC](https://github.com/s-andrews/FastQC) and adaptors trimmed with [TrimeGalore](https://github.com/FelixKrueger/TrimGalore). QC was repeated on trimmed reads to verify.

**Read Alignment & Read Counts** : BWA was used to index the [D39 reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014365.2/) and reads were then mapped back to it. [SAMtools](https://doi.org/10.1093/gigascience/giab008) were used to convert SAM->BAM and [featureCounts](https://doi.org/10.1093/bioinformatics/btt656) used to generate the gene-level count matrix. 

**Normalization** : To make gene expression values comparable across samples, normalization was done in three steps:
* `deseq_normalization`: [DeSeq2](https://doi.org/10.1093/gigascience/giab008) was used to remove library size effects and make counts comparable across samples -> later on used for differential expression analyses & PCA.
* `edger_normalization`: [edgeR](https://doi.org/10.1093/bioinformatics/btp616) was used to adjust for library compoistion effects using TMM (Trimmed Mean of M-values) -> later on used for exploratory analyses & clustering.
* `rpkm_calcluation`: normalize counts for gene length and sequencing depth -> compare expression levels across genes (not just across samples).

**Differntial Expression** : Differntially expressed genes between conditions (`data/conditions.tsv`). [DeSeq2](https://doi.org/10.1093/gigascience/giab008) was used with input files: normalized counts & conditions.tsv file.

The piepeline followed can be found in: [RNASeq](https://github.com/adamd3/BactSeq)
