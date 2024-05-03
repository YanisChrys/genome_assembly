############  PRIMARY PREPROCESSING MODULE  ############
############  IF READS ARE NOT CIRCULAR  ############

# Programs and Descriptions in current file:
# - CCS (Circular Consensus Sequencing):
#   - Program: ccs
#   - Generates highly accurate single-molecule sequences from PacBio subreads.

# - ACTC (Alignment-based Compression of Tandem Copies):
#   - Compresses tandem copies in subreads to improve downstream processing.

# - DeepConsensus:
#   - Uses machine learning to improve consensus accuracy in PacBio subreads.

# - Adapter Filtering or Renaming:
#   - Description: Filters out adapter sequences or renames files based on user preference.

# Conda Environments:
# - deepconsensus_ccs.yaml: Conda environment for running the CCS program.
# - actc.yaml: Conda environment for running the ACTC program.
# - adapterfilt.yaml: Conda environment for running the adapter filtering script.



# use pacbio subreads to create hifireads
# break up the reads into several pieces and
# run ccs in multiple parallel chunks to speed up the process
# afterwards, use deepconsensus to correct errors
# deepconsensus model can ben downloaded from 
# https://console.cloud.google.com/storage/browser/brain\
#-genomics-public/research/deepconsensus/models/v1.2/model_checkpoint?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))
rule ccs:
    input:
        config["subreads"]
    output:
        bamo=temp("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam"),
        index=temp("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam.pbi"),
        metrics="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}_metrics.json.gz",
        rprt="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}_report.txt"
    params:
        loglevel=config["ccs"]["loglevel"],
        chunk="{chunknumber}/%s" % config["ccs"]["chunks"]
    threads:
        config["ccs"]["cores"]
    log:
        "RESULTS/LOG/PREPROCESSING.CCS_PACBIO." + PREFIX + "_{chunknumber}.ccs.log"
    conda:
        "../envs/deepconsensus_ccs.yaml"
    shell: """
        ccs {input} {output.bamo} \
        --min-rq=0.95 \
        --num-threads {threads} \
        --chunk {params.chunk} \
        --log-level {params.loglevel} --log-file {log} --report-file {output.rprt} --metrics-json {output.metrics}
    """

# align subreads to CCS reads.
rule actc:
    input:
        ori_subs=config["subreads"],
        ccs_subs="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam"
    output:
        touch("RESULTS/PREPROCESSING/CCS_PACBIO/ACTC/" + PREFIX + "_{chunknumber}.subreads_to_ccs_actc.bam")
    threads:
        config["ccs"]["cores"]
    conda:
        "../envs/actc.yaml"
    priority: 1
    shell: """
        actc -j {threads}  \
        {input.ori_subs} \
        {input.ccs_subs} \
        {output}
    """

#  create new corrected reads in FASTQ format
rule deepconsensus:
    input:
        actc_subs="RESULTS/PREPROCESSING/CCS_PACBIO/ACTC/" + PREFIX + "_{chunknumber}.subreads_to_ccs_actc.bam",
        ccs_subs="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam"
    output:
        temp("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.fastq")
    params:
        model=config["modules"]["deepconsensus"]["model_checkpoint"],
        deepconsensus_run_method=config["modules"]["deepconsensus"]["run_method"]
    threads:
        min(workflow.cores,15)
    shell: """
        {params.deepconsensus_run_method}

        deepconsensus run \
        --subreads_to_ccs={input.actc_subs}  \
        --ccs_bam={input.ccs_subs} \
        --checkpoint={params.model} \
        --output={output} \
        --batch_zmws 50 \
        --cpus {threads}
    """

rule merge_deepconsensus_files:
    input:
        expand("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.fastq",chunknumber=CHUNK_NMB)
    output:
        "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/" + PREFIX + "_merged_ccs.fastq"
    shell: """
        cat {input} > {output}
    """

# you should be able to do adapter filtering with hifi from the run or from the config file
# optionally create a blast db of your adapters
rule filter_adapters:
    input:
        lambda wildcards: "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/" + PREFIX + "_merged_ccs.fastq" if DO_CCS else config["hifi"]
    output:
        "DATA/" + PREFIX + ".filtered.fastq.gz"
    conda:
        "../envs/adapterfilt.yaml"
    priority: 1
    params:
        my_prefix=lambda wildcards, input: os.path.basename(input[0]).replace('.fq.gz', '').replace('.fastq.gz', '')
    threads:
        workflow.cores
    shell: """
            #makeblastdb -in $PWD/DATA/HiFiAdapterFilt/DB/* -dbtype nucl
            export PATH=$PATH:$PWD/DATA/HiFiAdapterFilt/DB

            cp -n {input} DATA/
            cp -n utils/pbadapterfilt.sh DATA/
            cd DATA

            bash pbadapterfilt.sh -p {params.my_prefix} -l 45 -t {threads}

            mv -n {params.my_prefix}.filt.fastq.gz ../{output}
            rm -f $(basename {input})
    """
    