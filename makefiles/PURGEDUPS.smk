############ Purge erroneous duplications ############

# split fasta file where Ns occur

rule SPLIT_AT_Ns:
    input:
        "RESULTS/GENOME_ASSEMBLY/{prefix}_primary.fasta"
    output:
        "RESULTS/GENOME_ASSEMBLY/{prefix}_primary_split.fasta"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.split_fa.log"
    shell: """
        split_fa {input} > {output}
    """

# xasm5=asm-to-ref mapping, for ~0.1 sequence divergence (assembly to self)
# -I=split index for every ~NUM input bases

rule MINIMAP2_GENOME_TO_SELF:
    input:
        "RESULTS/GENOME_ASSEMBLY/{prefix}_primary_split.fasta"
    output:
        "RESULTS/MINIMAP2/{prefix}_primary_genome_split.paf"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.minimap2a.log"
    params:
        i_split=config["INDEX_SPLIT"]
    shell: """
        minimap2 -I {params.i_split} -t {threads} -xasm5 -DP {input} {input} > {output}
    """

rule MINIMAP2_GENOME_TO_READS:
    input:
        splitf="RESULTS/GENOME_ASSEMBLY/{prefix}_primary_split.fasta",
        ccs_reads="DATA/{prefix}.fastq.gz"
    output:
        "RESULTS/MINIMAP2/{prefix}_primary.reads_split.paf"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.minimap2b.log"
    params:
        i_split=config["INDEX_SPLIT"]
    shell: """
        minimap2 -I {params.i_split} -t {threads} -x map-pb {input.splitf} {input.ccs_reads} > {output}
    """

#Calculate haploid/diploid coverage threshold and remove haplotype duplicates from assembly

rule PBCSTAT:
    input:
        "RESULTS/MINIMAP2/{prefix}_primary.reads_split.paf"
    output:
        basecov="RESULTS/PURGE_DUPS/{prefix}/PB.base.cov",
        stat="RESULTS/PURGE_DUPS/{prefix}/PB.stat" #, PB.cov.wig
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.pbcstat.log"
    params:
        prefix="RESULTS/PURGE_DUPS/{prefix}"
    shell: """
        pbcstat -O {params.prefix} {input}
    """

rule CALCUTS:
    input:
        "RESULTS/PURGE_DUPS/{prefix}/PB.stat"
    output:
        "RESULTS/PURGE_DUPS/{prefix}/cutoffs"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.calcuts.log"
    params:
        prefix="RESULTS/PURGE_DUPS/{prefix}"
    shell: """
        calcuts {input} > {output}
    """

rule PURGE_DUPS:
    input:
        basecov="RESULTS/PURGE_DUPS/{prefix}/PB.base.cov",
        cutoffs="RESULTS/PURGE_DUPS/{prefix}/cutoffs",
        genome_paf="RESULTS/MINIMAP2/{prefix}_primary_genome_split.paf"
    output:
        "RESULTS/PURGE_DUPS/{prefix}_dups.bed"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.purge_dups.log"
    shell: """
        purge_dups -2 -c {input.basecov} -T {input.cutoffs} {input.genome_paf} > {output}
    """

rule GET_SEQS:
    input:
        purge_duped="RESULTS/PURGE_DUPS/{prefix}_dups.bed",
        fasta="RESULTS/GENOME_ASSEMBLY/{prefix}_primary.fasta"
    output:
        "RESULTS/PURGE_DUPS/{prefix}_seqs_purged.hap.fa"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.get_seqs.log"
    params:
        prefix="RESULTS/PURGE_DUPS/{prefix}_seqs_purged"
    shell: """
        get_seqs -e -p {params.prefix} {input.purge_duped} {input.fasta}
    """

rule INDEX_SEQS:
    input:
        "RESULTS/PURGE_DUPS/{prefix}_seqs_purged.hap.fa"
    output:
        "RESULTS/PURGE_DUPS/{prefix}_seqs_purged.hap.fa.fai"
    shell: """
        samtools faidx {input}
    """
