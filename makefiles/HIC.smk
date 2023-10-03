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

rule MAP_HIC_TO_ASSEMBLY:
    input:
        fasta=lambda wildcards: get_input(wildcards, HIFIASM_FILES, HIFIASM_BASENAMES, 'basename'),
        hic_reads=lambda wildcards: get_input(wildcards, HIC_READ_FILES, HIC_BASENAMES,'hic')
    output:
        temp("RESULTS/HIC/RAW/{basename}_{hic}.bam")
    threads:
        workflow.cores
    params:
        chunk=config["CHUNKS"]
    conda:
        "../envs/bwa.yaml"
    shell: """
        bwa index {input.fasta}
        bwa mem -t{threads} -B8 {input.fasta} {input.hic_reads} | samtools view -b - > {output}
    """

# Filter for 5’-most alignment for each read, scripts are salsa 2 scripts

rule FILTER_5PRIME_ALIGNMENT:
    input:
        rawhic1="RESULTS/HIC/RAW/{basename}_{hic}.bam"
    output:
        filtrawhic1=temp("RESULTS/HIC/FILTERED/{basename}_{hic}_filtered.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.{basename}_{hic}.filter_r1.log"
    params:
        chunk=config["CHUNKS"]
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools view -h  {input.rawhic1} | perl /share/scientific_bin/salsa2/2.3/bin/filter_five_end.pl | \
        samtools view -@ {threads} -b - > {output.filtrawhic1}
    """

# Combine read 1 and 2 maps with mapq>10 (!)
# use the script from the arima assembly to combine the hic reads
# uses the .fai index of our assembly to index the hic data
# perl script:
# <read 1 bam> <read 2 bam> <path to samtools> <minimum map quality filter>
rule COMBINE_HICR1_AND_HICR2:
    input:
        hicfiles=expand("RESULTS/HIC/FILTERED/{basename}_{hic}_filtered.bam", hic=HIC_BASENAMES,basename=HIFIASM_BASENAMES),
        faidxfile="RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa.fai"
    output:
        temp("RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads.bam")
    threads:
        workflow.cores
    conda:
        "../envs/combine_perls.yaml"
    envmodules:
        config["SALSA2_MODULE"]
    shell: """
        perl /share/scientific_bin/salsa2/2.3/bin/two_read_bam_combiner.pl {input.hicfiles} samtools 10 | \
        samtools view -@ {threads} -b -t {input.faidxfile} - > {output}
    """

# --RGLB \ -LB    null    Read-Group library
# --RGPL \ -PL    null    Read-Group platform (e.g. ILLUMINA, SOLID)
# --RGPU \ -PU    null    Read-Group platform unit (eg. run barcode)
# --RGSM \ -SM    null    Read-Group sample name
rule ADD_READ_GROUP:
    input:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads.bam"
    output:
        temp("RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_withRG.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC_{basename}.combine.log"
    params:
        prefix=PREFIX,
        label=SAMPLE_LABEL
    envmodules:
        config["GATK_MODULE"]
    shell: """ 
        gatk AddOrReplaceReadGroups \
        -I {input} -O {output} \
        -ID {params.prefix} -LB {params.label} -SM {params.label} -PL ILLUMINA -PU none
    """

rule SORT_HIC:
    input:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_withRG.bam"
    output:
        sorted=temp("RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_sorted.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.COMBINED.DEDUPED_{basename}.samtoolssort1.log"
    params:
        chunk=config["CHUNKS"]
    conda:
        "../envs/samtools.yaml"
    shell: """
        mkdir RESULTS/HIC/COMBINED/DEDUPED/TEMP/
        samtools sort -@ {threads} -m 1G -O bam -o {output.sorted} {input}
    """


rule REMOVE_DUPLICATE_MAPPINGS:
    input:
        sorted="RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_sorted.bam"
    output:
        deduped=temp("RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_deduped.bam"),
        metrics="RESULTS/HIC/COMBINED/DEDUPED/{basename}_picard_hic_deduped_metrics.txt"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC_{basename}.picdedup.log"
    envmodules:
        config["GATK_MODULE"]
    shell: """
        #-XX:MaxPermSize=1g -XX:+CMSClassUnloadingEnabled  
        gatk MarkDuplicates -XX:ParallelGCThreads={threads} -Xmx2g --REMOVE_DUPLICATES true \
        -I {input.sorted} -O {output.deduped} -M {output.metrics} \
        --ASSUME_SORT_ORDER "coordinate" --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1024 --SORTING_COLLECTION_SIZE_RATIO 0.1 --MAX_RECORDS_IN_RAM 250000
    """

rule RE_SORT_HIC:
    input:
        deduped="RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_deduped.bam"
    output:
        resorted=touch("RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_deduped_resorted.bam")
    threads:
        workflow.cores
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools sort -@ {threads} -n -m 1G -O bam -o {output.resorted} {input.deduped}
    """

############ YAHS ############

# input: contig assembly indexed fasta +  BAM/BED/BIN file with the alignment results of Hi-C reads to the contigs
# -l:  minimum contig length included for scaffolding.
# -q: minimum read mapping quality 

# -e {params.enzymes} add this with an if statement when enzymes is given
# YaHS can only handle up to 45,000 contigs. Consider excluding short contigs from scaffolding (with -l option) if the contig number exceeds this limit.
rule YAHS:
    input:
        hic_bam="RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_deduped_resorted.bam",
        fasta="RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa",
        fai="RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa.fai"
    output:
        "RESULTS/HIC/YAHS/{basename}_yahs_scaffolds_final.fa"
    threads:
        workflow.cores
    params:
        dir="RESULTS/HIC/YAHS/{basename}_yahs"
    conda:
        "../envs/yahs.yaml"
    shell: """
        yahs {input.fasta} {input.hic_bam} -o {params.dir} 
    """
