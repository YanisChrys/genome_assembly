

import re
import os
configfile: "config_stats.yaml"

# separate basename from extension and then only keep basename
ASSEMBLY_FILES = config["genomes"]
ASSEMBLY_BASENAMES = [os.path.splitext(os.path.basename(f))[0] for f in ASSEMBLY_FILES]

# split your assemplies at the alternate/primary so you're left with the basal 
# prefix of your assembly which be used as the identifier of the two haplotypes
MERQURY_BASENAMES = list(set(
    [os.path.basename(f).rsplit('_primary', 1)[0].rsplit('_alternate', 1)[0] for f in ASSEMBLY_FILES]
))

PREFIX=re.sub('[^a-zA-Z0-9 \_]', '', config["file_prefix"])

ASSEMBLY_RUN=config["run"]

def get_input(wildcards, files, basenames, wildcard_name):
    # For each pair, it checks if the basename matches the wildcard. If it does, it yields the file path.
    # within its parentheses. If no item satisfies the condition, a StopIteration exception is raised.

    return next(
        (f for f, b in zip(files, basenames) if b == getattr(wildcards, wildcard_name))
    )


rule all:
    input:
        expand("RESULTS/BUSCO/{basename}/run_{db}/missing_busco_list.tsv", db = config["busco"]["lineage"], basename = ASSEMBLY_BASENAMES),
        expand("RESULTS/BUSCO/{basename}/short_summary.specific.{db}.{basename}.txt", db = config["busco"]["lineage"], basename = ASSEMBLY_BASENAMES),
        expand("RESULTS/STATISTICS_POL/GFASTATS/{basename}_gfastats_report.txt", basename = ASSEMBLY_BASENAMES),
        expand("RESULTS/STATISTICS_POL/QUAST/{basename}/quast.log", basename = ASSEMBLY_BASENAMES),
        expand("RESULTS/STATISTICS/MERYL/{basename}_"+ASSEMBLY_RUN+".qv", basename = MERQURY_BASENAMES),
        expand("RESULTS/STATISTICS_POL/SEQKIT/{basename}.assemblyStat.txt", basename = ASSEMBLY_BASENAMES)




# STATS #

rule meryl:
    input:
        "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"
    output:
        directory("RESULTS/STATISTICS/MERYL/" + PREFIX + ".meryl")
    threads:
        min(workflow.cores,15)
    priority: 1
    params:
        wd="RESULTS/STATISTICS/MERYL/",
        out=PREFIX + ".meryl",
        kmers=config["meryl"]["kmer"]
    conda:
        "../envs/meryl.yaml"
    shell: """
        cd {params.wd}
        meryl count k={params.kmers} ../../../{input} threads={threads} memory=200g output {params.out}
    """

############  GENOME SIZE ESTIMATE WITH GENESCOPE.FK  ############

rule meryl_hist:
    input:
        "RESULTS/STATISTICS/MERYL/" + PREFIX + ".meryl"
    output:
        "RESULTS/STATISTICS/MERYL/" + PREFIX + ".hist"
    conda:
        "../envs/meryl.yaml"
    threads:
        1
    shell:
        "meryl histogram {input} > {output}"


rule seqkit_assembly_stats:
    input:
        lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/SEQKIT/{basename}.assemblyStat.txt"
    priority: 1
    conda:
        "envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """

rule merqury:
    input:
        meryldb="RESULTS/STATISTICS/MERYL/" + PREFIX + ".meryl",
        asm1=lambda wildcards: get_input_hifiasm(wildcards, ASSEMBLY_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[0],
        asm2=lambda wildcards: get_input_hifiasm(wildcards, ASSEMBLY_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[1]
    output:
        "RESULTS/STATISTICS/MERYL/{basename}_"+ASSEMBLY_RUN+".qv"
    threads:
        min(workflow.cores,15)
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        read_prefix=PREFIX,
        merqury_prefix="{basename}_hifiasm"
    shell: """
        cd RESULTS/STATISTICS/MERYL/
        export OMP_NUM_THREADS={threads}
        merqury.sh {input.meryldb} ../../../{input.asm1} ../../../{input.asm2} {params.merqury_prefix}
    """

rule gfastats:
    input:
        fasta=lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/GFASTATS/{basename}_gfastats_report.txt"
    conda:
        "envs/gfastas.yaml"
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
        "envs/quast.yaml"
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
        "envs/busco.yaml"
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
