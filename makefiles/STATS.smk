############  INITIALIZE PRIMARY STATISTICAL MODULE AND TABULAR ARRAY ############

# Pipeline Programs and purpose:
# - SEQKIT: Used for initial read statistics from the CCS reads in rules SEQKIT_CCS_STATS, SEQKIT_ASSEMBLY_STATS, and SEQKIT_YAHS_STATS.
# - FastK: Used for k-mer analysis in rule FASTK_TABLE.
# - Histex: Used for simplifying FastK hist table in rule HISTEX.
# - GeneScopeFK: Used for generating plots describing genome properties in rule GENESCOPE_FK.
# - KatGC: Used for generating 3D heat map or contour map of the frequency of a k-mer versus its' GC content in rule KATGC.
# - MerquryFK: Used for evaluating the completeness of the assembly in rule MerquryFK.
# - PloidyPlot: Used for generating ploidy plots in rule PLOIDYPLOT.
# - GAAS: Used for generating assembly statistics in rule gaas.
# - QUAST: Used for evaluating assembly quality in rule quast.
# - GT SEQSTATS: Used for generating sequence statistics in rule gtseqstat.
# - GFASTATS: Used for generating GFA statistics in rule gfastats.


rule SEQKIT_HIFI_STATS:
    input:
        lambda wildcards: "DATA/" + PREFIX + ".fastq.gz" if SUBREADS else config["HIFI"]
    output:
        "RESULTS/STATISTICS/SEQKIT/" + PREFIX + ".ccs.readStats.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """


rule FASTQ_HIC_STATS:
    input:
        "RESULTS/HIC/YAHS/" + PREFIX + "_yahs_scaffolds_final.fa"
    output:
        "RESULTS/STATISTICS/SEQKIT/" + PREFIX + "._yahs_scaffolds_final_assemblyStat.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
"""


# FastK
# kmers are k-sized string chunks - here it is DNA
# FastK splits the DNA into kmers and counts how many timers it sees each one
# need to provide absolute path of the module's locations
# See the README.mk for the prerequisite steps to running the programs from Fastk and Merqury packages
rule FASTK:
    input:
        lambda wildcards: "DATA/" + PREFIX + ".fastq.gz" if SUBREADS else config["HIFI"]
    output:
        ktab=touch("RESULTS/STATISTICS/FASTK/" + PREFIX + ".ktab"),
        hist=touch("RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist")
    threads:
        workflow.cores
    priority: 1
    params:
        out="RESULTS/STATISTICS/FASTK/" + PREFIX,
        kmers=config["KMER"]
    envmodules:
        config["FASTK_MODULE"]
    shell: """
        FastK -v -t1 -T{threads} -k{params.kmers} -N{params.out} {input}
    """

############  GENOME SIZE ESTIMATE WITH GENESCOPE.FK  ############

# Histex histograg=input for Genescope.FK
# Simplify FastK hist table
# h=frequencies to be displayed
# NOTE: consider exploring alterantive h values
rule HISTEX:
    input:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist"
    output:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist.txt"
    priority: 1
    log:
        "RESULTS/LOG/STATISTICS.FASTK." + PREFIX + ".histex.log"
    envmodules:
        config["FASTK_MODULE"]
    shell: """
        Histex -G {input} -h1000 > {output}
    """


# Genescope.FK
#  plots describing genome properties such as genome size, heterozygosity, and repetitiveness
# it will fail if there isn't enough data
rule GENESCOPE_FK:
    input:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist.txt"
    output:
        progress=touch("RESULTS/STATISTICS/FASTK/genescopeFK/" + PREFIX + "_progress.txt"),
        summary="RESULTS/STATISTICS/FASTK/genescopeFK/" + PREFIX + "_summary.txt"
    threads:
        workflow.cores
    priority: 1
    log:
        "RESULTS/LOG/STATISTICS.FASTK." + PREFIX + ".GeneScopeFK.log"
    params:
        inputprefix=PREFIX,
        outfolder="RESULTS/STATISTICS/FASTK/genescopeFK",
        kmers=config["KMER"]
    envmodules:
        config["FASTK_MODULE"],
        config["GENOMESCOPE_MODULE"]
    shell: """
        GeneScopeFK.R -i {input} -o {params.outfolder} -k {params.kmers} --ploidy 2 --num_rounds 4 -n {params.inputprefix}
    """

############  READ ANALYSIS  ############
##########  MERQURY.FK MODULE  ##########

# KatGC
# 3D heat map or contour map of the frequency of a k-mer versus its' GC content
# input=FastK ktab file (auto-detected)
# needs to move to directory where it will find the files
rule KATGC:
    input:
        hist="RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist",
        ktab="RESULTS/STATISTICS/FASTK/" + PREFIX + ".ktab"
    output:
        png1="RESULTS/STATISTICS/FASTK/" + PREFIX + ".st.png",
        png2="RESULTS/STATISTICS/FASTK/" + PREFIX + ".fi.png",
        png3="RESULTS/STATISTICS/FASTK/" + PREFIX + ".ln.png"
    threads:
        workflow.cores
    priority: 1
    conda:
        "../envs/merqury.yaml"
    envmodules:
        config["FASTK_MODULE"]
    params:
        my_prefix=PREFIX
    shell: """
        cd RESULTS/STATISTICS/FASTK/
        KatGC -T{threads} {params.my_prefix} {params.my_prefix}
    """

rule MerquryFK:
    input:
        hist="RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist",
        ktab="RESULTS/STATISTICS/FASTK/" + PREFIX + ".ktab",
        asm=lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS/FASTK/{basename}.completeness.stats"
    threads:
        workflow.cores
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        merqury_dir=config['MERQURY_DIR']
    envmodules:
        config["FASTK_MODULE"]
    shell: """
        export PATH=$PATH:{params.merqury_dir} 
        export PATH=$PATH:{params.merqury_dir}/*
        cd RESULTS/STATISTICS/FASTK/
        MerquryFK -T{threads} {wildcards.basename} ../../../{input.asm} {wildcards.basename}
    """

rule PLOIDYPLOT:
    input:
        hist="RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist",
        ktab="RESULTS/STATISTICS/FASTK/" + PREFIX + ".ktab"
    output:
        "RESULTS/STATISTICS/FASTK/PLOIDYPLOT/" + PREFIX + ".png"
    threads:
        min(workflow.cores,5)
    priority: 1
    params:
        merqury_dir=config['MERQURY_DIR'],
        my_prefix=PREFIX
    envmodules:
        config["FASTK_MODULE"]
    conda:
        "../envs/merqury.yaml"
    shell: """
        export PATH=$PATH:{params.merqury_dir} 
        export PATH=$PATH:{params.merqury_dir}/*
        cd RESULTS/STATISTICS/FASTK/
        PloidyPlot -T{threads} -vk -pdf -oPLOIDYPLOT/{params.my_prefix} {params.my_prefix}
    """

# ASSEMBLY STATS

rule SEQKIT_ASSEMBLY_STATS:
    input:
        lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS/SEQKIT/{basename}.assemblyStat.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """

rule gaas:
    input:
        lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS/GAAS/{basename}_gaasstats.txt"
    priority: 1
    conda:
        "../envs/gaas.yaml"
    shell: """
        gaas_fasta_statistics.pl -f {input} > {output}
    """

# QUAST

rule quast:
    input:
        contigs = lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS/{basename}_hifiasm/quast/quast.log"
    threads: 
        min(workflow.cores,20)
    priority: 1
    params:
        "RESULTS/STATISTICS/{basename}_hifiasm/quast"
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

# GT SEQSTATS


rule extractgenomeSize:
    input:
        genescope_output="RESULTS/STATISTICS/FASTK/genescopeFK/" + PREFIX + "_summary.txt"
    output:
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_genome_size.txt"
    run:
        genome_size = extract_max_genome_haploid_length(input.genescope_output)
        with open(output.genome_size_file, 'w') as f:
            f.write(str(genome_size))

rule gtseqstat:
    input:
        assembly= lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename'),
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_genome_size.txt"
    output:
        stats="RESULTS/STATISTICS/GTSEQSTATS/{basename}_gtseqstat.txt"
    priority: 1
    params:
        genome_size=lambda wildcards, input: open(input.genome_size_file).read().strip()
    conda:
        "../envs/genometools.yaml"
    shell: """
        # module load genometools/1.6.2
        gt seqstat -contigs -genome {params.genome_size} \
        {input.assembly} > {output.stats}
    """

# GFASTATS

rule gfastats:
    input:
        fasta= lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename'),
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_genome_size.txt"
    output:
        "RESULTS/STATISTICS/GFASTATS/{basename}_gfastats_report.txt"
    conda:
        "../envs/gfastas.yaml"
    threads:
        min(workflow.cores,10)
    priority: 1
    params:
        genome_size=lambda wildcards, input: open(input.genome_size_file).read().strip()
    shell: """
        gfastats \
        --input-sequence {input.fasta} \
        {params.genome_size} \
        --tabular \
        --threads {threads} \
        --out-bubbles \
        --nstar-report \
        --seq-report \
        --out-size scaffolds > {output}
    """
