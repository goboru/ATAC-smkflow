# Snakefile

# By Goboru, November 2025
# Pipeline to run the complete pipeline to analyze ATAC data or do individual steps 


# Use: conda activate ATACpipeline; snakemake -s ATAC_pipeline.smk --cores 8

# Configurable paths
configfile: "config.yaml"  # optional

# CONFIGURATION
dir_raw = config.get("raw_dir", "raw_data")
dir_out = config.get("out_dir", "results")

# Gather FASTQ files
SAMPLES, = glob_wildcards(f"{dir_raw}/{{sample}}.fastq.gz")
UNIQ_SAMPLES, = glob_wildcards(f"{dir_raw}/{{uniq_sample}}_r1.fastq.gz")

# This rule sets when the pipeline is finished
rule all:
    input: 
        expand(f"{dir_out}/mito/{{uniq_sample}}.noMT.bam", uniq_sample=UNIQ_SAMPLES)
        #expand(f"{dir_out}/temp_trimming/{{sample}}.trimmed.fastq.gz", sample=SAMPLES)


# Rule 1: FastQC before trimming
rule fastqc_untrimmed:
    input:
        f"{dir_raw}/{{sample}}.fastq.gz"
    output:
        html=f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.html",
        zip=f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.zip"
    log:
        f"{dir_out}/logs/qc_untrimmed/{{sample}}.log"
    shell:
        "fastqc {input} --outdir {dir_out}/qc_untrimmed &> {log}"

# Rule 1.2: multiqc of FastQC's before trimming
rule multiqc_untrimmed:
    input:
        expand(f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.zip", sample=SAMPLES)
    output:
        html=f"{dir_out}/qc_untrimmed/multiqc_report.html"
    shell:
        "multiqc {dir_out}/qc_untrimmed --force -o {dir_out}/qc_untrimmed"




# Rule 2: Trimming
rule trimming:
    input:
        r1 = f"{dir_raw}/{{uniq_sample}}_r1.fastq.gz",
        r2 = f"{dir_raw}/{{uniq_sample}}_r2.fastq.gz"
    params:
        mode = config["trim_mode"]
    output:
        r1 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r1.trimmed.fastq.gz",
        r2 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r2.trimmed.fastq.gz",
        json = f"{dir_out}/temp_trimming/{{uniq_sample}}.fastp.json",
        html = f"{dir_out}/temp_trimming/{{uniq_sample}}.fastp.html"
    log:
        f"{dir_out}/logs/trimming/{{uniq_sample}}.log"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            {params.mode} \
            -j {output.json} -h {output.html} \
            &> {log}
        """


# Rule 3: FastQC after trimming
rule fastqc_trimmed:
    input:
        f"{dir_out}/temp_trimming/{{sample}}.trimmed.fastq.gz"
    output:
        html=f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.html",
        zip=f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.zip"
    log:
        f"{dir_out}/logs/qc_trimmed/{{sample}}.log"
    shell:
        "fastqc {input} --outdir {dir_out}/qc_trimmed &> {log}"

# Rule 3.2: multiqc of FastQC's after trimming
rule multiqc_trimmed:
    input:
        expand(f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.zip", sample=SAMPLES)
    output:
         html=f"{dir_out}/qc_trimmed/multiqc_report.html"
    shell:
        "multiqc {dir_out}/qc_trimmed --force -o {dir_out}/qc_trimmed"


# Rule 4: Alignment with bowtie2

rule bowtie2_align:
    input:
        r1 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r1.trimmed.fastq.gz",
        r2 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r2.trimmed.fastq.gz"
    output:
        bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
    log:
        bowtie = f"{dir_out}/logs/bowtie2/{{uniq_sample}}.log",
        samtools = f"{dir_out}/logs/samtools_sort/{{uniq_sample}}.log"
    params:
        index = config["bowtie2_index"]
    threads: config["bowtie2_threads"]
    shell:
        """
        (bowtie2 --very-sensitive -I 25 -X 700 \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} \
            2> {log.bowtie}) \
        | samtools sort \
            -@ 2 \
            -o {output.bam} \
            - \
            2> {log.samtools}
        """
        




# section 5: Post-alingments QCs

#Rule 5.1: idxstats reports
rule idxstats:
    input:
        bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
    output:
        report = f"{dir_out}/idx_report/{{uniq_sample}}.idxstats.txt"
    log:
        f"{dir_out}/logs/idx_report/{{uniq_sample}}.idxstats.log"
    threads: 1
    shell:
        """
        # Make sure the BAM is indexed
        samtools index {input.bam} 2>> {log}

        # Generate idxstats report
        samtools idxstats {input.bam} > {output.report} 2>> {log}
        """

# Rule 5.2: remove mitochondrial reads
rule remove_mito:
    input:
        bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
    output:
        bam = f"{dir_out}/mito/{{uniq_sample}}.noMT.bam"
    log:
        f"{dir_out}/logs/mito/{{uniq_sample}}.remove_mito.log"
    threads: 4
    shell:
        r"""
        samtools view -h {input.bam} \
            | awk '$3 != "chrM" && $3 != "MT" || $1 ~ /^@/' \
            | samtools view -b -o {output.bam} \
            &> {log}
        samtools index {output.bam}
        """

# Rule 5.3: Remove ENCODE blacklist regions
rule remove_blacklist:
    input:
        bam = f"{dir_out}/filtered/{{uniq_sample}}.noMT.bam"
    output:
        bam = f"{dir_out}/filtered/{{uniq_sample}}.noMT.noBlacklist.bam"
    log:
        f"{dir_out}/logs/filter/{{uniq_sample}}.remove_blacklist.log"
    params:
        blacklist = config["blacklist_bed"]
    threads: 4
    shell:
        r"""
        bedtools intersect \
            -v \
            -abam {input.bam} \
            -b {params.blacklist} \
            > {output.bam} 2> {log}

        samtools index {output.bam}
        """
# Remove non canonical reads?

# Rule 5.4: Mark duplicates, remove dupliccates and low quality reads
rule mark_duplicates:
    input:
        bam = f"{dir_out}/filtered/{{uniq_sample}}.noMT.noBlacklist.bam"
    output:
        bam  = f"{dir_out}/filtered/{{uniq_sample}}.final.dedup.bam",
        metrics = f"{dir_out}/filtered/{{uniq_sample}}.final.dedup.metrics.txt"
    log:
        f"{dir_out}/logs/filter/{{uniq_sample}}.dedup.log"
    threads: 4
    shell:
        r"""
        picard MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT \
            2> {log}

        samtools index {output.bam}
        """

# Rule : BAM samtools flagstat


# Rule 6: Peak calling with MACS2

# Rule 7: Quality controls after peak calling


# Rule 8: Differential analysis
# Rule 8.1: Sensitivity analysis. Sex, batch, age?