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
rule CCS:
    input:
        config["SUBREADS"]
    output:
        bamo=temp("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam"),
        index=temp("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam.pbi"),
        metrics="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}_metrics.json.gz",
        rprt="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}_report.txt"
    params:
        loglevel=config["LOGLEVEL"],
        chunk="{chunknumber}/%s" % config["CHUNKS"]
    threads:
        config["CORES"]
    log:
        "RESULTS/LOG/PREPROCESSING.CCS_PACBIO." + PREFIX + "_{chunknumber}.ccs.log"
    conda:
        "../envs/deepconsensus_ccs.yaml"
    shell: """
        ccs {input} {output.bamo} \
        --min-rq=0.88 \
        --num-threads {threads} \
        --chunk {params.chunk} \
        --log-level {params.loglevel} --log-file {log} --report-file {output.rprt} --metrics-json {output.metrics}
    """

rule ACTC:
    input:
        ori_subs=config["SUBREADS"],
        ccs_subs="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam"
    output:
        touch("RESULTS/PREPROCESSING/CCS_PACBIO/ACTC/" + PREFIX + "_{chunknumber}.subreads_to_ccs_actc.bam")
    threads:
        config["CORES"]
    conda:
        "../envs/actc.yaml"
    shell: """
        actc -j {threads}  \
        {input.ori_subs} \
        {input.ccs_subs} \
        {output}
    """

rule DEEPCONSENSUS:
    input:
        actc_subs="RESULTS/PREPROCESSING/CCS_PACBIO/ACTC/" + PREFIX + "_{chunknumber}.subreads_to_ccs_actc.bam",
        ccs_subs="RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam"
    output:
        temp("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.fastq")
    params:
        model=config["DEEPCONSENSUS_MODEL_CHECKPOINT"],
        deepconsensus_run_method=config["DEEPCONSENSUS_RUN_METHOD"]
    threads: 
        15
    shell: """
        {params.deepconsensus_run_method}

        deepconsensus run \
        --subreads_to_ccs={input.actc_subs}  \
        --ccs_bam={input.ccs_subs} \
        --checkpoint={params.model} \
        --output={output} \
        --cpus {threads}
    """

rule MERGE_DEEPCONSENSUS_FILES:
    input:
        expand("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.fastq",chunknumber=CHUNK_NMB, prefix=PREFIX)
    output:
        "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/" + PREFIX + ".merged.ccs.fastq"
    shell: """
        cat {input} > {output}
    """

# you should be able to do adapter filtering with hifi from the run or from the config file
rule FILTER_ADAPTERS_OR_RENAME:
    input:
        lambda wildcards: "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/" + PREFIX + ".merged.ccs.fastq" if DO_CCS else config["HIFI"]
    output:
        "DATA/" + PREFIX + ".fastq.gz"
        # filtfq=temp("DATA/" + PREFIX + ".merged.ccs.filt.fastq.gz") if DO_ADAPT_FILT else "DATA/" + PREFIX + ".fastq.gz"
    conda:
        "../envs/adapterfilt.yaml"
    priority: 1
    params:
        my_basename=PREFIX + ".merged.ccs",
        my_prefix=PREFIX,
        adapt_filt=str(DO_ADAPT_FILT).lower()  # Convert the boolean to a lowercase string ("true" or "false")
    shell: """
        if [ "{params.adapt_filt}" = "true" ]; then
            #makeblastdb -in $PWD/DATA/HiFiAdapterFilt/DB/* -dbtype nucl
            export PATH=$PATH:$PWD/DATA/HiFiAdapterFilt/DB
            cd RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/
            bash ../../../../utils/pbadapterfilt.sh -p {params.my_basename} -o ../../../../DATA/
            cp ../../../../DATA/{params.my_prefix}.merged.ccs.filt.fastq.gz ../../../../{output}
            rm ../../../../DATA/{params.my_prefix}.merged.ccs.filt.fastq.gz
        else
            mv -n {input} {output}
        fi
    """
