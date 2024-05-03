# BUSCO Score

# -m: genome mode
# -i: input
# -o: output folder
# -l: lineage/database - define it in config file
# write a script for busco
#write a rule that downloads the environment and make it a local one

#"RESULTS/BUSCO/" + PREFIX + "/run_{db}/missing_busco_list.tsv"

rule busco:
    input:
        lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        summary="RESULTS/BUSCO/{basename}/short_summary.specific.{db}.{basename}.txt",
        missing="RESULTS/BUSCO/{basename}/run_{db}/missing_busco_list.tsv"
    threads:
        min(workflow.cores,20)
    priority: 1
    params:
        dataset_dir=config["busco"]["dataset_folder"],
        out_dir="RESULTS/BUSCO/",
        run_name="{basename}",
        lineage=config["busco"]["lineage"],
        plot_wd="RESULTS/BUSCO/{basename}/"
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
    """

rule busco_plots:
    input: 
        expand("RESULTS/BUSCO/{basename}/short_summary.specific.{db}.{basename}.txt",basename=ASSEMBLY_BASENAMES, db=config["busco"]["lineage"])
    output: 
        "RESULTS/BUSCO/busco_figure.png"
    params:
        plot_wd="RESULTS/BUSCO/"
    conda:
        "../envs/busco.yaml"
    shell: """
        cp -n {input} {params.plot_wd}

        generate_plot.py -q -wd {params.plot_wd}
    """

