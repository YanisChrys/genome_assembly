############  HIC ASSEMBLY  ############

# Takes a an assembly fasta file "" and maps

############  HIC SCAFFOLDING  ############

# Map HiC reads to contigs with Salsa2
# the - symbol after samtools is to tell samtools to take the input from the pipe
rule MAP_HIC_TO_ASSEMBLY:
    input:
        fasta=expand("RESULTS/PURGE_DUPS/{prefix}_seqs_purged.hap.fa",prefix=FILE_PREFIX),
        hic1="DATA/{hic}.fastq.gz"
    output:
        temp(expand("RESULTS/HIC/RAW/{prefix}_{hic}.bam",prefix=FILE_PREFIX,allow_missing=True))
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.map_bwa_{hic}.log"
    params:
        chunk=config["CHUNKS"]
    shell: """
        bwa index {input.fasta}
        bwa mem -t24 -B8 {input.fasta} {input.hic1} | samtools view -b - > {output}
    """

# Filter for 5’-most alignment for each read, scripts are salsa 2 scripts

rule FILTER_5PRIME_ALIGNMENT:
    input:
        rawhic1=expand("RESULTS/HIC/RAW/{prefix}_{hic}.bam",prefix=FILE_PREFIX,allow_missing=True)
    output:
        filtrawhic1=temp(expand("RESULTS/HIC/FILTERED/{prefix}_{hic}_filtered.bam",prefix=FILE_PREFIX,allow_missing=True))
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.{hic}.filter_r1.log"
    params:
        chunk=config["CHUNKS"]
    shell: """
        samtools view -h  {input.rawhic1} | perl /share/scientific_bin/salsa2/2.3/bin/filter_five_end.pl | \
        samtools view -@ {threads} -b - > {output.filtrawhic1}
    """

# Combine read 1 and 2 maps with mapq>10 (!)
# use the script from the arima assembly to combine the hic reads
# uses the .fai index of our assembly to index the hic data
rule COMBINE_HICR1_AND_HICR2:
    input:
        hicfiles=expand("RESULTS/HIC/FILTERED/{prefix}_{hic}_filtered.bam",hic=["hic1","hic2"],prefix=FILE_PREFIX),
        faidxfile=expand("RESULTS/PURGE_DUPS/{prefix}_seqs_purged.hap.fa.fai",prefix=FILE_PREFIX)
    output:
        temp("RESULTS/HIC/COMBINED/mapped_hic_reads.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.combine.log"
    params:
        chunk=config["CHUNKS"]
    # envmodules:
    #     "/share/scientific_bin/salsa2/2.3"
    shell: """
        module load salsa2/2.3
        perl /share/scientific_bin/salsa2/2.3/bin/two_read_bam_combiner.pl {input.hicfiles} samtools 10 | \
        samtools view -@ {threads} -b -t {input.faidxfile} - > {output}
    """

# replace dedup.sh from salsa2: this is better because it's more snakemake-like
# --RGLB \ -LB	null	Read-Group library
# --RGPL \ -PL	null	Read-Group platform (e.g. ILLUMINA, SOLID)
# --RGPU \ -PU	null	Read-Group platform unit (eg. run barcode)
# --RGSM \ -SM	null	Read-Group sample name
rule ADD_READ_GROUP:
    input:
        "RESULTS/HIC/COMBINED/mapped_hic_reads.bam"
    output:
        temp("RESULTS/HIC/COMBINED/mapped_hic_reads_withRG.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.combine.log"
    params:
        chunk=config["CHUNKS"],
        prefix=FILE_PREFIX,
        label=SAMPLE_LABEL
    # envmodules:
    #     "/share/scientific_bin/salsa2/2.3"
    shell: """
        picard AddOrReplaceReadGroups \
        -I {input} -O {output} \
        -ID {params.prefix} -LB {params.label} -SM {params.label} -PL ILLUMINA -PU none
    """

# Filter duplicate mappings and sort by name
# dedup.sh missing from salsa2 module on cluster
rule SORT_HIC:
    input:
        "RESULTS/HIC/COMBINED/mapped_hic_reads_withRG.bam"
    output:
        sorted=temp("RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_sorted.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.COMBINED.DEDUPED.samtoolssort1.log"
    params:
        chunk=config["CHUNKS"],
        tempf="RESULTS/HIC/COMBINED/DEDUPED/TEMP/mapped_hic_reads_withRG_sort.bam.tmp"
    shell: """
        mkdir RESULTS/HIC/COMBINED/DEDUPED/TEMP/
        samtools sort -@ {threads} -T {params.tempf} -m 1G -O bam -o {output.sorted} {input}
    """


rule REMOVE_DUPLICATE_MAPPINGS:
    input:
        sorted="RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_sorted.bam"
    output:
        deduped=temp("RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_deduped.bam"),
        metrics="RESULTS/HIC/COMBINED/DEDUPED/picard_hic_deduped_metrics.txt"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.picdedup.log"
    params:
        tempf="RESULTS/HIC/COMBINED/DEDUPED/TEMP/"
    shell: """
        #-XX:MaxPermSize=1g -XX:+CMSClassUnloadingEnabled
        picard MarkDuplicates -XX:ParallelGCThreads={threads} -Xmx2g -Djava.io.tmpdir={params.tempf} --REMOVE_DUPLICATES true \
        -I {input.sorted} -O {output.deduped} -M {output.metrics} --TMP_DIR {params.tempf} \
        --ASSUME_SORT_ORDER "coordinate" --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1024 --SORTING_COLLECTION_SIZE_RATIO 0.1 --MAX_RECORDS_IN_RAM 250000
    """

rule RE_SORT_HIC:
    input:
        deduped="RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_deduped.bam"
    output:
        resorted=touch("RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_deduped_resorted.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.samtoolssort2.log"
    params:
        chunk=config["CHUNKS"],
        tempf="RESULTS/HIC/COMBINED/DEDUPED/TEMP/"
    shell: """
        samtools sort -@ {threads} -n -T {params.tempf} -m 1G -O bam -o {output.resorted} {input.deduped}
    """

# Convert to bed file - salsa 2 desired input

rule BAM2BED_HIC:
    input:
        "RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_deduped_resorted.bam"
    output:
        temp("RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_deduped_resorted.bed")
    log:
        "RESULTS/LOG/HIC.bam2bed.log"
    params:
        chunk=config["CHUNKS"]
    shell: """
        bedtools bamtobed -i {input} > {output}
    """

############ SCAFFOLDING ############

# Run salsa2 scaffolding
rule RUN_SALSA2:
    input:
        hic_bed="RESULTS/HIC/COMBINED/DEDUPED/mapped_hic_reads_withRG_deduped_resorted.bed",
        fasta=expand("RESULTS/PURGE_DUPS/{prefix}_seqs_purged.hap.fa",prefix=FILE_PREFIX),
        fai=expand("RESULTS/PURGE_DUPS/{prefix}_seqs_purged.hap.fa.fai",prefix=FILE_PREFIX)
    output:
        agp="RESULTS/HIC/SALSA2/scaffolds_FINAL.agp",
        fastout=protected("RESULTS/HIC/SALSA2/scaffolds_FINAL.fasta")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/HIC.salsa2.log"
    params:
        chunk=config["CHUNKS"],
        enzymes=config["ENZYMES"],
        dir="RESULTS/HIC/SALSA2/"
    conda:
        "envs/salsa2.yaml"
    shell: """
        #source /share/scientific_bin/anaconda3/2022.05/etc/profile.d/conda.sh
        #conda activate salsa2
        python2 /share/scientific_bin/salsa2/2.3/bin/run_pipeline.py --assembly {input.fasta} --length {input.fai} \
        --enzyme {params.enzymes} --bed {input.hic_bed} --output {params.dir} \
        --clean yes --prnt yes --iter 50
    """


rule RUN_BUSCO_AGAIN:
    input:
        "RESULTS/HIC/SALSA2/scaffolds_FINAL.fasta"
    output:
        # prefix="RESULTS/BUSCO/{prefix}_",
        # short_json="RESULTS/BUSCO/{prefix}_short_summary.json",
        # short_txt="RESULTS/BUSCO/{prefix}_short_summary.txt",
        # full_table="RESULTS/BUSCO/{prefix}_full_table.tsv",
        # touch("RESULTS/BUSCO/{prefix}_busco_missing.tsv")
        # directory("RESULTS/BUSCO_SCAFFOLDED/{prefix}/")
        "RESULTS/BUSCO_SCAFFOLDED/{prefix}/run_{db}/missing_busco_list.tsv"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}_{db}.busco.scaf.log"
    params:
        dataset_dir="busco_downloads",
        out_dir="RESULTS/BUSCO_SCAFFOLDED/",
        run_name=FILE_PREFIX,
        chunk=config["CHUNKS"],
        lineage=config["LINEAGE"]
    shell: """
        busco -m genome \
        -c {threads} \
        -i {input} \
        -o {params.run_name}
        --download_path {params.dataset_dir} \
        --out_path {params.out_dir} \
        -l {params.lineage} \
        --offline
    """