############  INITIALIZE PRIMARY STATISTICAL MODULE AND TABULAR ARRAY ############

# initial read statistics from the CCS reads
# Do all (-a) statistics on the reads (quartiles of seq length, sum_gap, N50, etc)


rule SEQKIT_CCS_STATS:
    input:
        "DATA/{prefix}.fastq.gz"
    output:
        "RESULTS/GENOME_ASSEMBLY/PREPROCESSING/STATISTICS/{prefix}.ccs.readStats.txt"
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/PREPROCESSING.STATISTICS.{prefix}.seqkit.log"
    shell: """
        seqkit stats -a {input} > {output}
    """

# FastK
# kmers are k-sized string chunks - here it is DNA
# FastK splits the DNA into kmers and counts how many timers it sees each one
# need to provide absolute path of the module's locations
# See the README.mk for the prerequisite steps to running the programs from Fastk and Merqury packages
rule FASTK_TABLE:
    input:
        "DATA/{prefix}.fastq.gz"
    output:
        ktab=touch("RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.ktab"),
        hist=touch("RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.hist")
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/PREPROCESSING.STATISTICS.{prefix}.fastk.log"
    params:
        out="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}",
        temp="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/tmp",
        chunk=config["CHUNKS"]
    # envmodules:
    #     "fastk/current"
    shell: """
        module load fastk/current
        FastK -v -t1 -T{threads} -k40 -N{params.out} -P{params.temp} {input}
    """
############  GENOME SIZE ESTIMATE WITH GENESCOPE.FK  ############

# Histex histograg=input for Genescope.FK
# Simplify FastK hist table
# h=frequencies to be displayed
# NOTE: consider exploring alterantive h values
rule HISTEX:
    input:
        "RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.hist"
    output:
        "RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.hist.txt"
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/PREPROCESSING.STATISTICS.{prefix}.histex.log"
    # envmodules:
    #     "fastk/current"
    shell: """
        module load fastk/current
        Histex -G {input} -h1000 > {output}
    """

# Genescope.FK
#  plots describing genome properties such as genome size, heterozygosity, and repetitiveness
# it will fail if there isn't enough data
rule GENESCOPE_FK:
    input:
        "RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.hist.txt"
    output:
        progress=touch("RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/genescopeFK/{prefix}_progress.txt")
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/PREPROCESSING.STATISTICS.{prefix}.GeneScopeFK.log"
    params:
        inputprefix="{prefix}",
        outfolder="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/genescopeFK",
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"],
        kmersize=config["KMER"]
    shell: """
        module load genescopefk/current
        module load fastk/current
        GeneScopeFK.R -i {input} -o {params.outfolder} -k {params.kmersize} -n {params.inputprefix}
    """

############  READ ANALYSIS  ############
##########  MERQURY.FK MODULE  ##########

# module path: /share/scientific_bin/merquryfk/current

# KatGC
# 3D heat map or contour map of the frequency of a k-mer versus its' GC content
# input=FastK ktab file (auto-detected)
# needs to move to directory where it will find the files
rule KATGC:
    input:
        hist="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.hist",
        ktab="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.ktab"
    output:
        png1="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.st.png",
        png2="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.fi.png",
        png3="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.ln.png"
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/PREPROCESSING.STATISTICS.{prefix}.KatGC.log"
    envmodules:
        "merquryfk/current",
        "fastk/current"
    conda:
        "envs/merqury.yaml"
    shell: """
        # module load merquryfk/current
        # module load fastk/current
        # source /share/scientific_bin/anaconda3/2022.05/etc/profile.d/conda.sh
        # conda activate merqury
        cd RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/
        KatGC -T{threads} {wildcards.prefix} {wildcards.prefix}
    """

# PloidyPlot
# useful in case the ploidy of your sample is unknown
# input=FastK ktab file (auto-detected)
# verbose and keep the table it creates
# plot in pdf format
# Needs program Logex
rule PLOIDYPLOT:
    input:
        hist="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.hist",
        ktab="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.ktab"
    output:
        touch("RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/{prefix}.pdf")
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/PREPROCESSING.STATISTICS.{prefix}.PloidyPlot.log"
    params:
        inp=lambda wildcards, output: output[0][:-4],
        out="RESULTS/GENOME_ASSEMBLY/PREPROCESSING/FASTKMERS/",
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"]
    envmodules:
        "merquryfk/current",
        "fastk/current"
    shell: """
        PloidyPlot -T{threads} -v -pdf -o{params.out} {params.inp}
    """