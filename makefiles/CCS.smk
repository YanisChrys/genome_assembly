############  PRIMARY PREPROCESSING MODULE  ############
############  IF READS ARE NOT CIRCULAR  ############

# run ccs to create circular hifi reads

# use pacbio subreads to create hifireads
# break up the reads into several piecies and
# run ccs in multiple parallel chunks to speed up the process
rule RUN_CCS:
    input:
        "DATA/{prefix}.subreads.bam"
    output:
        bamo="RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}.ccs.bam",
        index="RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}.ccs.bam.pbi",
        metrics="RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}_metrics.json.gz",
        rprt="RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}_report.txt"
    threads:
        config["CORES"]
    params:
        loglevel=config["LOGLEVEL"],
        chunk="{chunknumber}/%s" % config["CHUNKS"],
    log:
        "RESULTS/LOG/PREPROCESSING.CCS_PACBIO.{prefix}_{chunknumber}.ccs.log"
    shell: """
        ccs {input} {output.bamo} --num-threads {threads} --chunk {params.chunk} \
        --log-level {params.loglevel} --log-file {log} --report-file {output.rprt} --metrics-json {output.metrics}
    """

rule MERGE_CCS:
    input:
        expand("RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}.ccs.bam",chunknumber=CHUNK_NMB, prefix=FILE_PREFIX)
    output:
        touch("RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/{prefix}.merged.ccs.bam")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/PREPROCESSING.CCS_PACBIO.{prefix}.merge.log"
    shell: """
        samtools merge -c -p -f -@{threads} {output} {input}
    """

rule INDEX_MERGED_CCS:
    input:
        "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/{prefix}.merged.ccs.bam"
    output:
        touch("RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/{prefix}.merged.ccs.bam.pbi")
    log:
        "RESULTS/LOG/PREPROCESSING.CCS_PACBIO.{prefix}.pbindex.log"
    shell: """
        pbindex {input}
    """

rule CONVERT_TO_FASTQ:
    input:
        bam="RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/{prefix}.merged.ccs.bam",
        pbi="RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/{prefix}.merged.ccs.bam.pbi"
    output:
        "DATA/{prefix}.fastq.gz"
    params:
        "DATA/{prefix}"
    log:
        "RESULTS/LOG/PREPROCESSING/CCS_PACBIO/{prefix}.bam2fastq.log"
    shell: """
        bam2fastq -o {params} {input.bam}
    """
