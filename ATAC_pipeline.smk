# Snakefile

# By Goboru, November 2025
# Pipeline to run the complete pipeline to analyze ATAC data or do individual steps 

# Use: snakemake -s ATAC_pipeline.smk --cores 8 --use-conda

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
        expand(f"{dir_out}/aligned/{{uniq_sample}}_align.bam", uniq_sample=UNIQ_SAMPLES)
        #f"{dir_out}/qc_trimmed/multiqc_report.html"
        #expand(f"{dir_out}/temp_trimming/{{sample}}.trimmed.fastq.gz", sample=SAMPLES)
        #expand(f"{dir_out}/aligned/{{sample}}_align.bam", sample=UNIQ_SAMPLES)

#Rule 0: Create directory for saving the logs
rule make_log_dir:
    output:
        directory("logs")
    shell:
        "mkdir -p logs"

# Rule 1: FastQC before trimming
rule fastqc_untrimmed:
    input:
        f"{dir_raw}/{{sample}}.fastq.gz"
    output:
        html=f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.html",
        zip=f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.zip"
    log:
        f"logs/qc_untrimmed/{{sample}}.log"
    shell:
        "fastqc {input} --outdir {dir_out}/qc_untrimmed &> {log}"

# Rule 1.2: multiqc of FastQC's before trimming
rule multiqc_untrimmed:
    input:
        expand(f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.zip", sample=SAMPLES)
    output:
        f"{dir_out}/qc_untrimmed/multiqc_report.html"
    shell:
        "multiqc {dir_out}/qc_untrimmed --force --outdir {dir_out}/qc_untrimmed"




# Rule 2: Trimming
rule trimming:
    input:
        r1 = f"{dir_raw}/{{uniq_sample}}_r1.fastq.gz",
        r2 = f"{dir_raw}/{{uniq_sample}}_r2.fastq.gz"
    output:
        r1 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r1.trimmed.fastq.gz",
        r2 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r2.trimmed.fastq.gz",
        json = f"{dir_out}/temp_trimming/{{uniq_sample}}.fastp.json",
        html = f"{dir_out}/temp_trimming/{{uniq_sample}}.fastp.html"
    log:
        f"logs/trimming/{{uniq_sample}}.log"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            -j {output.json} -h {output.html} \
            &> {log}
        """


# Rule 3: FastQC before trimming
rule fastqc_trimmed:
    input:
        f"{dir_out}/temp_trimming/{{sample}}.trimmed.fastq.gz"
    output:
        html=f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.html",
        zip=f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.zip"
    log:
        f"logs/qc_trimmed/{{sample}}.log"
    shell:
        "fastqc {input} --outdir {dir_out}/qc_trimmed &> {log}"

# Rule 3.2: multiqc of FastQC's after trimming
rule multiqc_trimmed:
    input:
        expand(f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.zip", sample=SAMPLES)
    output:
        f"{dir_out}/qc_trimmed/multiqc_report.html"
    shell:
        "multiqc {dir_out}/qc_trimmed --force --outdir {dir_out}/qc_trimmed"



# Rule 4: Alignment with bowtie2
# rule bowtie2_align:
#     input:
#         r1 = f"{dir_raw}/{{uniq_sample}}_r1.fastq.gz",
#         r2 = f"{dir_raw}/{{uniq_sample}}_r2.fastq.gz"
#     output:
#         bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
#     log:
#         err = f"logs/bowtie2/{{uniq_sample}}.err",
#         out = f"logs/bowtie2/{{uniq_sample}}.out"
#     params:
#         index = config["bowtie2_index"]
#     threads: 6
#     shell:
#         """
#         bowtie2 --very-sensitive -I 25 -X 700 -k 10 \
#             -x {params.index} \
#             -1 {input.r1} -2 {input.r2} \
#             -p {threads} \
#             2> {log.err} | samtools sort -o {output.bam} - \
#             > {log.out}
#         """
        
rule bowtie2_sam:
    input:
        r1 = f"{dir_raw}/{{uniq_sample}}_r1.fastq.gz",
        r2 = f"{dir_raw}/{{uniq_sample}}_r2.fastq.gz"
    output:
        sam = f"{dir_out}/aligned/{{uniq_sample}}.sam"
    log:
        f"logs/bowtie2/{{uniq_sample}}.log"
    params:
        index = config["bowtie2_index"]
    threads: 6
    shell:
        """
        bowtie2 --very-sensitive -I 25 -X 700 -k 10 \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} \
            1> {output.sam} \
            2> {log}
        """


# Rule 4.1: Samtools after bowtie2
    #This should have been the same step, I could not make it work the two in one

rule samtools_sort:
    input:
        sam = f"{dir_out}/aligned/{{uniq_sample}}.sam"
    output:
        bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
    log:
        f"logs/samtools_afterbw2/{{uniq_sample}}.log"
    threads: 4
    shell:
        """
        samtools sort -@ {threads} \
            -o {output.bam} \
            {input.sam} \
            &> {log}
        """

# Rule 5: Post-alingments QCs

# Rule 5.1: remove mitochondrial reads
# #Generate the idxstats report
# samtools idxstats <sample>_sorted.bam > <sample>_sorted.idxstats

# #Check the number of reads mapped to the mitochondria (chrM)
# grep "chrM" <sample>_sorted.idxstats
# #Generate the flagstat report
# samtools flagstat <sample>_sorted.bam > <sample>_sorted.flagstat

# #check the total number of aligned fragments
# head <sample>_sorted.flagstat
# #Remove reads aligned to the mitochondria
# samtools view -h <sample>_sorted.bam | grep -v chrM | samtools sort -O bam -o <sample>.rmChrM.bam -T .

# Rule 5.2: Mark duplicates, remove dupliccates and low quality reads

# Rule 5.3: Remove ENCODE blacklist regions


# Rule 6: Peak calling with MACS2

# Rule 7: Quality controls after peak calling

# Rule 8: Differential analysis
