# BUSCO Score

# -m: genome mode
# -i: input
# -o: output folder
# -l: lineage/database - define it in config file
# write a script for busco
#write a rule that downloads the environment and make it a local one
rule RUN_BUSCO:
    input:
        "RESULTS/GENOME_ASSEMBLY/{prefix}_primary.fasta"
    output:
        # prefix="RESULTS/BUSCO/{prefix}_",
        # short_json="RESULTS/BUSCO/{prefix}_short_summary.json",
        # short_txt="RESULTS/BUSCO/{prefix}_short_summary.txt",
        # full_table="RESULTS/BUSCO/{prefix}_full_table.tsv",
        # touch("RESULTS/BUSCO/{prefix}_busco_missing.tsv")
        # directory("RESULTS/BUSCO/{prefix}/")
        "RESULTS/BUSCO/{prefix}/run_{db}/missing_busco_list.tsv"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}_{db}.busco.log"
    params:
        dataset_dir="busco_downloads",
        out_dir="RESULTS/BUSCO/",
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
