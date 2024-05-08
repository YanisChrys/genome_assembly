# STATS #

rule fastk:
    input:
        config["reads"]
    output:
        ktab=touch("RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".ktab"),
        hist=touch("RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist")
    threads:
        workflow.cores
    priority: 1
    params:
        out="RESULTS/STATISTICS_POL/FASTK/" + PREFIX,
        kmers=config["modules"]["fastk"]["kmer"]
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        FastK -v -t1 -T{threads} -k{params.kmers} -N{params.out} {input}
    """

rule histex:
    input:
        "RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist"
    output:
        "RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist.txt"
    priority: 1
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        Histex -G {input} > {output}
    """

rule MerquryFK:
    input:
        hist="RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist",
        ktab="RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".ktab",
        asm=lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/FASTK/{basename}.completeness.stats"
    threads:
        workflow.cores
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        merqury_dir=config["modules"]["merqury"]["path"],
        my_prefix=PREFIX
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        export PATH=$PATH:{params.merqury_dir} 
        export PATH=$PATH:{params.merqury_dir}/*
        cd RESULTS/STATISTICS_POL/FASTK/
        MerquryFK -T{threads} {params.my_prefix} ../../../{input.asm} {wildcards.basename}
    """


rule seqkit_assembly_stats:
    input:
        lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/SEQKIT/{basename}.assemblyStat.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """

rule gfastats:
    input:
        fasta=lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/GFASTATS/{basename}_gfastats_report.txt"
    conda:
        "../envs/gfastas.yaml"
    threads:
        min(workflow.cores,10)
    priority: 1
    params:
        genome_size=config["genome_size"]
    shell: """
        gfastats \
        {input.fasta} \
        {params.genome_size} \
        --tabular \
        --threads {threads} \
        --nstar-report \
        --seq-report \
        --out-size scaffolds > {output}
    """


rule quast:
    input:
        contigs = lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/QUAST/{basename}/quast.log"
    threads: 
        min(workflow.cores,20)
    priority: 1
    params:
        "RESULTS/STATISTICS_POL/QUAST/{basename}"
    conda:
        "../envs/quast.yaml"
    shell: """
        quast.py \
        --eukaryote \
        --large \
        -t {threads} \
        --min-contig 300 \
        -o {params} \
        {input.contigs}
    """

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
