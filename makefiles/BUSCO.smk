# BUSCO Score

# -m: genome mode
# -i: input
# -o: output folder
# -l: lineage/database - define it in config file
# write a script for busco
#write a rule that downloads the environment and make it a local one

#"RESULTS/BUSCO/" + PREFIX + "/run_{db}/missing_busco_list.tsv"

rule BUSCO:
    input:
        lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/BUSCO/{basename}/{basename}/run_{db}/missing_busco_list.tsv"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{basename}_{db}.busco.log"
    params:
        dataset_dir=config["BUSCO_DATASET_FOLDER"],
        out_dir="RESULTS/BUSCO/",
        run_name=PREFIX,
        lineage=config["LINEAGE"],
        plot_wd="RESULTS/BUSCO/{basename}/{basename}/"
    conda:
        "../envs/busco.yaml"
    shell: """
        busco -m genome -f \
        -c {threads} \
        -i {input} \
        -o {params.run_name} \
        --download_path {params.dataset_dir} \
        --out_path {params.out_dir} \
        -l {params.lineage} \
        --offline

        generate_plot.py -q -wd {params.plot_wd}
    """

