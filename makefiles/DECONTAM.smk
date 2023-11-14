
rule KRAKEN2:
    input:
        lambda wildcards: "DATA/" + PREFIX + ".fastq.gz" if SUBREADS else config["HIFI"]
    output:
        classification="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classification.txt",
        report="DATA/DECONTAMINATED/" + PREFIX + "_kraken.report.csv",
        decontaminated="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz",
        contamination="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq.gz"
    threads:
        workflow.cores
    conda:
        "../envs/kraken2.yaml"
    params:
        krakendb=config["KRAKEN_DB"],
        conf=config["KRAKEN_CONFIDENCE"]
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
        
        module load pigz/2.7

        pigz -p {threads} {params.decontaminated} && \
        pigz -p {threads} {params.contamination}
    """

