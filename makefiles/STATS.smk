############  INITIALIZE PRIMARY STATISTICAL MODULE AND TABULAR ARRAY ############

# initial read statistics from the CCS reads
# Do all (-a) statistics on the reads (quartiles of seq length, sum_gap, N50, etc)


rule SEQKIT_CCS_STATS:
    input:
        "DATA/{prefix}.fastq.gz"
    output:
        "RESULTS/STATISTICS/SEQKIT/{prefix}.ccs.readStats.txt"
    log:
        "RESULTS/LOG/STATISTICS.FASTK.{prefix}.seqkit.log"
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
        ktab=touch("RESULTS/STATISTICS/FASTK/{prefix}.ktab"),
        hist=touch("RESULTS/STATISTICS/FASTK/{prefix}.hist")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/STATISTICS.FASTK.{prefix}.fastk.log"
    params:
        out="RESULTS/STATISTICS/FASTK/{prefix}",
        temp="RESULTS/STATISTICS/FASTK/tmp",
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
        "RESULTS/STATISTICS/FASTK/{prefix}.hist"
    output:
        "RESULTS/STATISTICS/FASTK/{prefix}.hist.txt"
    log:
        "RESULTS/LOG/STATISTICS.FASTK.{prefix}.histex.log"
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
        "RESULTS/STATISTICS/FASTK/{prefix}.hist.txt"
    output:
        progress=touch("RESULTS/STATISTICS/FASTK/genescopeFK/{prefix}_progress.txt")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/STATISTICS.FASTK.{prefix}.GeneScopeFK.log"
    params:
        inputprefix="{prefix}",
        outfolder="RESULTS/STATISTICS/FASTK/genescopeFK",
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
        hist="RESULTS/STATISTICS/FASTK/{prefix}.hist",
        ktab="RESULTS/STATISTICS/FASTK/{prefix}.ktab"
    output:
        png1="RESULTS/STATISTICS/FASTK/{prefix}.st.png",
        png2="RESULTS/STATISTICS/FASTK/{prefix}.fi.png",
        png3="RESULTS/STATISTICS/FASTK/{prefix}.ln.png"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/STATISTICS.FASTK.{prefix}.KatGC.log"
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
        cd RESULTS/STATISTICS/FASTK/
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
        hist="RESULTS/STATISTICS/FASTK/{prefix}.hist",
        ktab="RESULTS/STATISTICS/FASTK/{prefix}.ktab"
    output:
        touch("RESULTS/STATISTICS/FASTK/{prefix}.pdf")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/STATISTICS.FASTK.{prefix}.PloidyPlot.log"
    params:
        inp=lambda wildcards, output: output[0][:-4],
        out="RESULTS/STATISTICS/FASTK/",
        chunk=config["CHUNKS"]
    envmodules:
        "merquryfk/current",
        "fastk/current"
    shell: """
        PloidyPlot -T{threads} -v -pdf -o{params.out} {params.inp}
    """
