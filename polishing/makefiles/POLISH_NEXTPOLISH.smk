

import re
import os
import glob
configfile: "config/config_nextpolish.yaml"

# get sample basename regardless of the extension (fa, fasta) and create sample wildcard
asm_path=config["asm_path"]

asm_path="/home/ychrysostomakis/EpiAus/asm_fasta/"
extensions = ["fa", "fasta","fa.gz", "fasta.gz"]
patterns = [os.path.join(asm_path, f"*.{ext}") for ext in extensions] # create list of possible targets
files = [file for pattern in patterns for file in glob.glob(pattern)] # now look up these combinations and find the true files
samples = {os.path.basename(file).rsplit('.', 1)[0] for file in files} # create the sample wildcard out of those true samples

samples=["EpiAus_202401_3SMRT_hifiasm_incpert_D10N150_l2_primary"]
# get the basename of the reads so fastk only runs once
PREFIX = re.sub(r'\.(fq\.gz|fastq\.gz|fq|fastq)$', '', os.path.basename(config["reads"]))

# characterist of your run to be added to outputs:
run=config["run"]

def get_input(wildcards, files, basenames, wildcard_name):
    return next(
        (f for f, b in zip(files, basenames) if b == getattr(wildcards, wildcard_name))
    )

FA_FILES = (expand("RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_" + run + "-1.fa", sample = samples) +
    expand("RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_" + run + "-2.fa", sample = samples) +
    expand("RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_" + run + "-3.fa", sample = samples))

BASENAMES = [os.path.basename(f).split('.')[0] for f in FA_FILES]

rule all:
    input:
    # polishing
        expand("RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_"+run+"-3.fa", sample = samples),
    # stats
        expand("RESULTS/BUSCO/{basename}/run_{db}/missing_busco_list.tsv", db = config["busco"]["lineage"], basename=BASENAMES),
        expand("RESULTS/BUSCO/{basename}/short_summary.specific.{db}.{basename}.txt", db = config["busco"]["lineage"], basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/GFASTATS/{basename}_gfastats_report.txt", basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/QUAST/{basename}/quast.log", basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/FASTK/{basename}.completeness.stats", basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/SEQKIT/{basename}.assemblyStat.txt", basename=BASENAMES)

# nextpolish
# create input config file for nextpolish and use for run
# tempdir should be temp: fix!
rule run_nextpolish:
    input:
        asm=lambda wildcards: next(os.path.join(asm_path, f"{wildcards.sample}.{ext}")
                               for ext in extensions
                               if glob.glob(os.path.join(asm_path, f"{wildcards.sample}.{ext}"))),
        pacbioreads=config["reads"]
    output:
        out1="RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_"+run+"-1.fa",
        out2="RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_"+run+"-2.fa",
        out3="RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_"+run+"-3.fa"
    threads: 
        workflow.cores
    conda:
        "../envs/nextpolish.yaml"
    envmodules:
        config["modules"]["nextpolish"]["bbmap"],
        config["modules"]["nextpolish"]["bwa"],
        config["modules"]["nextpolish"]["samtools"]
    params:
        genome_size=config["genome_size"],
        multithread_jobs=config["modules"]["nextpolish"]["multithread_jobs"],
        parallel_jobs=config["modules"]["nextpolish"]["parallel_jobs"],
        tempdir="RESULTS/ASM_POLISHING_NEXTPOLISH/{sample}_"+run+"-wkdir",
        exe_path=config["modules"]["nextpolish"]["path"]
    shell: """
        mkdir -p {params.tempdir}
        cp {input.asm} {params.tempdir}/asm.fa
        cd {params.tempdir}

        cat <<EOF > hifi.fofn 
{input.pacbioreads}
EOF

        ### NEXTPOLISH RUN 1
        cat <<EOF > run-1.cfg 
job_type = local
job_prefix = run-1
task = best
rewrite = yes
multithread_jobs = {params.multithread_jobs}
parallel_jobs = {params.parallel_jobs}
genome = ./asm.fa
genome_size = {params.genome_size}
workdir = ./01_rundir
polish_options = -p {{multithread_jobs}}
hifi_fofn = ./hifi.fofn 
hifi_minimap2_options = -x asm20 -t {{multithread_jobs}}
EOF
        
        {params.exe_path}/nextPolish ./run-1.cfg;
        cp ./01_rundir/genome.nextpolish.fasta ../../../{output.out1}


        ### NEXTPOLISH RUN 2
        cat <<EOF > ./run-2.cfg 
job_type = local
job_prefix = run-2
task = best
rewrite = yes
multithread_jobs = {params.multithread_jobs}
parallel_jobs = {params.parallel_jobs}
genome = ./01_rundir/genome.nextpolish.fasta
genome_size = {params.genome_size}
workdir = ./02_rundir
polish_options = -p {{multithread_jobs}}
hifi_fofn = ./hifi.fofn
hifi_minimap2_options = -x asm20 -t {{multithread_jobs}}
EOF
        
        {params.exe_path}/nextPolish ./run-2.cfg;
        cp ./02_rundir/genome.nextpolish.fasta ../../../{output.out2}
        
        ### NEXTPOLISH RUN 3
        cat  <<EOF > ./run-3.cfg
job_type = local
job_prefix = run-3
task = best
rewrite = yes
multithread_jobs = {params.multithread_jobs}
parallel_jobs = {params.parallel_jobs}
genome = ./02_rundir/genome.nextpolish.fasta
genome_size = {params.genome_size}
workdir = ./03_rundir
polish_options = -p {{multithread_jobs}}
hifi_fofn = ./hifi.fofn
hifi_minimap2_options = -x asm20 -t {{multithread_jobs}}
EOF
        
        {params.exe_path}/nextPolish ./run-3.cfg;
        cp ./03_rundir/genome.nextpolish.fasta ../../../{output.out3}
"""



# STATS #

rule fastk:
    input:
        reads=config["reads"],
        hic1=config["hic1"],
        hic2=config["hic2"]
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
        FastK -v -t1 -T{threads} -k{params.kmers} -N{params.out} {input.reads} {input.hic1} {input.hic2}
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
        asm=lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/FASTK/{basename}.completeness.stats"
    threads:
        workflow.cores
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        merqury_dir=config["modules"]["merqury"]["path"],
        my_prefix=PREFIX,
        outprefix="{basename}"
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        export PATH=$PATH:{params.merqury_dir} 
        export PATH=$PATH:{params.merqury_dir}/*
        cd RESULTS/STATISTICS_POL/FASTK/
        MerquryFK -T{threads} {params.my_prefix} ../../../{input.asm} {params.outprefix}
    """


rule seqkit_assembly_stats:
    input:
        lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
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
        fasta=lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
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
        contigs=lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
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
        lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
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

 