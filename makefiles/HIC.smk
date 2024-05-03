############  HIC ASSEMBLY  ############

# Programs and Descriptions in current file:
# - Map Hi-C reads to the assembly.
#   - Program: BWA and Samtools.

# - Filter for the 5’-most alignment for each read.
#   - Program: Samtools and a Perl script from SALSA2.

# - Combine read 1 and 2 maps with a minimum map quality of 10.
#   - Program: A Perl script from SALSA2 and Samtools.

# - Add read group information to the BAM file.
#   - Program: GATK.

# - Sort the Hi-C BAM file.
#   - Program: Samtools.

# - Filter out duplicate mappings from the BAM file.
#   - Program: GATK.

# - Re-sort the deduplicated Hi-C BAM file.
#   - Program: Samtools.

# - Scaffold the assembly using Hi-C data.
#   - Program: YaHS.
# -B8: penalty for mismatch 8 from 4 default
rule map_hic_to_assembly:
    input:
        fasta="RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta",
        hic_reads=lambda wildcards: get_input(wildcards, HIC_READ_FILES, HIC_BASENAMES,'hic')
    output:
        "RESULTS/HIC/RAW/{basename}_{hic}.bam"
    threads:
        min(workflow.cores,30)
    conda:
        "../envs/bwa.yaml"
    shell: """
        bwa-mem2 index {input.fasta}
        # increased mismatch penalty compared to default
        bwa-mem2 mem -t {threads} -B8 {input.fasta} {input.hic_reads} | \
        samtools view -@ {threads} -Sb -o {output} -
    """

# Filter for 5’-most alignment for each read, scripts are salsa 2 scripts

rule filter_5prime_alignment:
    input:
        rawhic1="RESULTS/HIC/RAW/{basename}_{hic}.bam"
    output:
        filtrawhic1="RESULTS/HIC/FILTERED/{basename}_{hic}_filtered.bam"
    threads:
        min(workflow.cores,15)
    log:
        "RESULTS/LOG/HIC.{basename}_{hic}.filter_r1.log"
    params:
        salsa_bin=config["modules"]["salsa2"]["bin_folder"]
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools view -h {input.rawhic1} | perl {params.salsa_bin}/filter_five_end.pl samtools 10 | \
        samtools view -@ {threads} -Sb - > {output.filtrawhic1}
    """

# Combine read 1 and 2 maps with mapq>10 (!)
# use the script from the arima assembly to combine the hic reads
# uses the .fai index of our assembly to index the hic data
# perl script:
# <read 1 bam> <read 2 bam> <path to samtools> <minimum map quality filter>
rule combine_hicr1_and_hicr2:
    input:
        hicfiles=[
            "RESULTS/HIC/FILTERED/{basename}_"+HIC_BASENAMES[0]+"_filtered.bam",
            "RESULTS/HIC/FILTERED/{basename}_"+HIC_BASENAMES[1]+"_filtered.bam"],
        faidxfile="RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta.fai"
    output:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads.bam"
    threads:
        min(workflow.cores,15)
    conda:
        "../envs/combine_perls.yaml"
    params:
        salsa_bin=config["modules"]["salsa2"]["bin_folder"]
    envmodules:
        config["modules"]["salsa2"]["path"]
    shell: """
        perl {params.salsa_bin}/two_read_bam_combiner.pl {input.hicfiles} samtools 30 | \
        samtools view -@ {threads} -Sb -t {input.faidxfile} - > {output}
    """

# --RGLB \ -LB    null    Read-Group library
# --RGPL \ -PL    null    Read-Group platform (e.g. ILLUMINA, SOLID)
# --RGPU \ -PU    null    Read-Group platform unit (eg. run barcode)
# --RGSM \ -SM    null    Read-Group sample name



rule sort_hic:
    input:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads.bam"
    output:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted.bam"
    threads:
        min(workflow.cores,5)
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools sort -@ {threads} -m 5G -O bam -o {output} {input}
    """

rule add_read_group:
    input:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted.bam"
    output:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG.bam"
    threads:
        min(workflow.cores,10)
    log:
        "RESULTS/LOG/HIC_{basename}.combine.log"
    params:
        prefix=PREFIX
    envmodules:
        config["modules"]["gatk"]["path"]
    shell: """ 
        gatk AddOrReplaceReadGroups --java-options "-Xmx5g" \
        -I {input} -O {output} \
        -ID {params.prefix} -LB lib_{params.prefix} -SM lib_{params.prefix} -PL ILLUMINA -PU none
    """

# SORTING_COLLECTION_SIZE_RATIO: help avoid running out of memory
# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: maximum number of open files
# MAX_RECORDS_IN_RAM: reduce amount of ram needed
rule remove_duplicate_mappings:
    input:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG.bam"
    output:
        deduped="RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG_deduped.bam",
        metrics="RESULTS/HIC/COMBINED/DEDUPED/{basename}_picard_hic_deduped_metrics.txt"
    threads:
        min(workflow.cores,10)
    log:
        "RESULTS/LOG/HIC_{basename}.picdedup.log"
    envmodules:
        config["modules"]["gatk"]["path"]
    shell: """ 
        gatk MarkDuplicates --java-options "-Xmx5g" --REMOVE_DUPLICATES true \
        -I {input} -O {output.deduped} -M {output.metrics} \
        --ASSUME_SORT_ORDER "coordinate" \
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1024 \
        --SORTING_COLLECTION_SIZE_RATIO 0.1 \
        --MAX_RECORDS_IN_RAM 250000
    """

# yahs input bam needs to be aligned by read name
rule re_sort_hic:
    input:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG_deduped.bam"
    output:
        resorted="RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG_deduped_resorted.bam",
        flagstat="RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_deduped_resorted_flagstat.tsv"
    threads:
        min(workflow.cores,15)
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools sort -@ {threads} -n -m 5G -O bam -o {output.resorted} {input}        
        samtools flagstat -O tsv {output.resorted} > {output.flagstat}
    """

############ YAHS ############

# input: contig assembly indexed fasta +  BAM/BED/BIN file with the alignment results of Hi-C reads to the contigs
# -l:  minimum contig length included for scaffolding.
# -q: minimum read mapping quality 

# -e {params.enzymes} add this with an if statement when enzymes is given
# YaHS can only handle up to 45,000 contigs. Consider excluding short contigs from scaffolding (with -l option) if the contig number exceeds this limit.
# yahs has an upper limit of 45000 contigs. Do not run if you have more
rule yahs:
    input:
        hic_bam="RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG_deduped_resorted.bam",
        fasta="RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta",
        fai="RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta.fai"
    output:
        "RESULTS/HIC/YAHS/{basename}_yahs_scaffolds_final.fa"
    threads:
        min(workflow.cores,20)
    priority: 1
    params:
        dir="RESULTS/HIC/YAHS/{basename}_yahs"
    conda:
        "../envs/yahs.yaml"
    shell: """
        yahs {input.fasta} {input.hic_bam} -o {params.dir} 
    """

# input for pretext
# hic stats, etc
rule map_scaf2hic:
    input:
        fasta=lambda wildcards: get_input(wildcards, SCAF_FILES, SCAF_BASENAMES, 'basename'),
        hic1=HIC1,
        hic2=HIC2
    output:
        "RESULTS/HIC/BAM4VIZ/{basename}.bam"
    threads:
        min(workflow.cores,60)
    conda:
        "../envs/bwa.yaml"
    params:
        rggroup = "'@RG\tID:{basename}\tSM:{basename}\tLB:{basename}\tPL:ILLUMINA'"
    shell: """
        bwa-mem2 index {input.fasta}
        bwa-mem2 mem -t{threads} -B8 -M -R {params.rggroup} {input.fasta} {input.hic1} {input.hic2} | \
        samtools view -@ {threads} -Sb - | samtools sort -@ {threads} -o {output} -
    """


