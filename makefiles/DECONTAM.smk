
rule KRAKEN2:
    input:
        lambda wildcards: "DATA/" + PREFIX + ".fastq.gz" if SUBREADS else config["HIFI"]
    output:
        classification="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classification.txt",
        report="DATA/DECONTAMINATED/" + PREFIX + "_kraken.report.csv",
        decontaminated="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq",
        contamination="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq"
    threads:
        workflow.cores
    conda:
        "../envs/kraken2.yaml"
    params:
        krakendb=config["KRAKEN_DB"]
        conf=config["KRAKEN_CONFIDENCE"]
    shell: """
        kraken2 --db {params.krakendb} \
            --threads {threads} \
            --use-names \
            --confidence  \
            --output {output.classification} \
            --report {output.report} \
            --unclassified-out {output.decontaminated} \
            --classified-out {output.contamination} \
            {input}
    """

