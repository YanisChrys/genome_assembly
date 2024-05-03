
rule kraken2:
    input:
        lambda wildcards: (
        "DATA/" + PREFIX + ".filtered.fastq.gz" if DO_ADAPT_FILT 
        else "RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/" + PREFIX + "_merged_ccs.fastq" if SUBREADS 
        else config["hifi"]
        )
    output:
        classification="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classification.txt",
        report="DATA/DECONTAMINATED/" + PREFIX + "_kraken.report.csv",
        decontaminated="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz",
        contamination="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq.gz"
    threads:
        min(workflow.cores,20)
    conda:
        "../envs/kraken2.yaml"
    params:
        krakendb=config["kraken"]["db"],
        conf=config["kraken"]["confidence"],
        decontaminated="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq",
        contamination="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq"
    shell: """
        kraken2 --db {params.krakendb} \
            --threads {threads} \
            --use-names \
            --confidence {params.conf} \
            --output {output.classification} \
            --report {output.report} \
            --unclassified-out {params.decontaminated} \
            --classified-out {params.contamination} \
            {input}

        pigz -p {threads} {params.decontaminated} && \
        pigz -p {threads} {params.contamination}
    """
