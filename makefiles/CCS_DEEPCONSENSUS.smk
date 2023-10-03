############  PRIMARY PREPROCESSING MODULE  ############
############  IF READS ARE NOT CIRCULAR  ############

# Remove reads with adapter sequences ("adapter contamination")
# from subreads


rule FILTER_ADAPTERS:
    input:
        "DATA/{prefix}.subreads.bam"
    output:
        "DATA/{prefix}.subreads_adapt_filt.bam"
    conda:
        "envs/HiFiAdapterFilt.yaml"
    shell: """
        bash ./utils/pbadapterfilt.sh 
    """

# run ccs to create circular hifi reads

# use pacbio subreads to create hifireads
# break up the reads into several pieces and
# run ccs in multiple parallel chunks to speed up the process
# afterwards, use deepconsensus to correct errors
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
        chunk="{chunknumber}/%s" % config["CHUNKS"]
    shell: """
        ccs {input} {output.bamo} \
        --min-rq=0.88 --num-threads {threads} --chunk {params.chunk} \
        --log-level {params.loglevel} --log-file {log} --report-file {output.rprt} --metrics-json {output.metrics}
    """


rule ACTC:
    input:
        ori_subs="DATA/{prefix}.subreads.bam",
        ccs_subs="RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}.ccs.bam"
    output:
        touch("RESULTS/PREPROCESSING/CCS_PACBIO/ACTC/{prefix}_{chunknumber}.subreads_to_ccs_actc.bam")
    threads:
        workflow.cores
    shell: """
        actc -j {threads}  \
        {input.ori_subs} \
        {input.ccs_subs} \
        {output}
    """

rule DEEPCONSENSUS:
    input:
        actc_subs="RESULTS/PREPROCESSING/CCS_PACBIO/ACTC/{prefix}_{chunknumber}.subreads_to_ccs_actc.bam",
        ccs_subs="RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}.ccs.bam"
    output:
        "RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}.ccs.fastq"
    envmodules:
        "deepconsensus/0.3.1"
    shell: """
        deepconsensus run \
        --subreads_to_ccs={input.actc_subs}  \
        --ccs_bam={input.ccs_subs} \
        --checkpoint=model/checkpoint \
        --output={output}
    """

rule MERGE_FQ:
    input:
        expand("RESULTS/PREPROCESSING/CCS_PACBIO/{prefix}_{chunknumber}.ccs.fastq",chunknumber=CHUNK_NMB, prefix=FILE_PREFIX)
    output:
        "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/{prefix}.merged.ccs.fastq"
    shell: """
        cat {input} > output
    """

rule COMPRESS_FQ:
    input:
        "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/{prefix}.merged.ccs.fastq"
    output:
        "DATA/{prefix}.fastq.gz"
    shell: """
        gzip -c {input} > {output}
    """
