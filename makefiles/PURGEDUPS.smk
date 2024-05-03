############ Purge erroneous duplications ############

# Programs and Descriptions in current file:
# Split the FASTA file wherever 'N' nucleotides occur 
# with a custom script (`split_fa`) from purge_dups package.

# - Map the assembly to itself to identify regions of simility with minimap2.

# - Map the split genome to HiFi reads with minimap2.

# - Calculate the base-level coverage for the mapped reads with `pbcstat`.

# - Determine coverage thresholds to differentiate between haploid and diploid regions with `calcuts`.

# - Purge Duplicates with `purge_dups`.

# - Extract sequences that are not marked as duplicates with get_seqs.

# - Index the purged sequences with `samtools faidx`


# split fasta file where Ns occur
# this way of generating the input has two functions:
# 1) it finds the correct file and 
# 2) it introduces the basename wildcard to the rule which can then be used downstream
rule split_at_ns:
    input:
        lambda wildcards: get_input(wildcards, HIFIASM_FILES, HIFIASM_BASENAMES, 'basename')
    output:
        temp("RESULTS/PURGE_DUPS/{basename}_split.fasta")
    threads:
        5
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        split_fa {input} > {output}
    """

# xasm5=asm-to-ref mapping, for ~0.1 sequence divergence (assembly to self)
# -I=split index for every ~NUM input bases

rule minimap2_genome_to_self:
    input:
        "RESULTS/PURGE_DUPS/{basename}_split.fasta"
    output:
        "RESULTS/MINIMAP2/{basename}_genome_split.paf"
    threads:
        min(workflow.cores,20)
    params:
        i_split=config["index_split"]
    conda:
        "../envs/minimap2.yaml"
    shell: """
        minimap2 -I {params.i_split} -t {threads} -x asm5 -DP {input} {input} > {output}
    """
#-x asm20
rule minimap2_genome_to_reads:
    input:
        splitf="RESULTS/PURGE_DUPS/{basename}_split.fasta",
        hifi_reads="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"
    output:
        "RESULTS/MINIMAP2/{basename}_reads_split.paf"
    threads:
        min(workflow.cores,20)
    params:
        i_split=config["index_split"]
    conda:
        "../envs/minimap2.yaml"
    shell: """
        minimap2 -I {params.i_split} -t {threads} -x map-hifi {input.splitf} {input.hifi_reads} > {output}
    """


#Calculate haploid/diploid coverage threshold and remove haplotype duplicates from assembly

rule pbcstat:
    input:
        "RESULTS/MINIMAP2/{basename}_reads_split.paf"
    output:
        basecov="RESULTS/PURGE_DUPS/{basename}/PB.base.cov",
        stat="RESULTS/PURGE_DUPS/{basename}/PB.stat" #, PB.cov.wig
    threads:
        5
    params:
        prefix="RESULTS/PURGE_DUPS/{basename}"
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        pbcstat -O {params.prefix} {input}
        #python3 scripts/hist_plot.py -c cutoff_file PB.stat PB.cov.png
    """

rule calcuts:
    input:
        "RESULTS/PURGE_DUPS/{basename}/PB.stat"
    output:
        "RESULTS/PURGE_DUPS/{basename}/cutoffs"
    threads:
        5
    params:
        prefix="RESULTS/PURGE_DUPS/{basename}"
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        calcuts {input} > {output}
    """
#  python3 hist_plot.py -c cutoffs PB.stat PB.cov.png
rule purge_dups:
    input:
        basecov="RESULTS/PURGE_DUPS/{basename}/PB.base.cov",
        cutoffs="RESULTS/PURGE_DUPS/{basename}/cutoffs",
        genome_paf="RESULTS/MINIMAP2/{basename}_genome_split.paf"
    output:
        temp("RESULTS/PURGE_DUPS/{basename}_dups.bed")
    threads:
        5
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        purge_dups -2 -c {input.basecov} -T {input.cutoffs} {input.genome_paf} > {output}
    """

rule get_seqs:
    input:
        purge_duped="RESULTS/PURGE_DUPS/{basename}_dups.bed",
        fasta="RESULTS/GENOME_ASSEMBLY/{basename}.fasta"
    output:
        "RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta"
    threads:
        5
    params:
        prefix="RESULTS/PURGE_DUPS/{basename}_seqs_purged",
        getseqsoutput="RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa"
    conda:
        "../envs/purge_dups.yaml"
    priority: 1
    shell: """
        get_seqs -e -p {params.prefix} {input.purge_duped} {input.fasta}
        mv {params.getseqsoutput} {output}
    """
#-e
rule index_seqs:
    input:
        "RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta"
    output:
        "RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta.fai"
    conda:
        "../envs/samtools.yaml"
    threads:
        5
    shell: """
        samtools faidx {input}
    """
