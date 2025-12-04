# ATAC-seq analysis with snakemake 

* Initial QC with fastqc and multiqc
* Trimming with fastp
* Post-trimming QC with fastqc and multiqc
* Alignment to hg38 with bowtie2
* Remove mitochondrial reads
* Remove ENCODE blacklist regions

Pending:
- Peak calling with macs2
- Remove mitochondrial reads
- Remove duplicate reads with picard
- Calculate FRiP score
- Calculate TSS enrichment and plot
- Generate bigWig tracks for visualization
- Diff accesibility analysis
