############ Purge erroneous duplications ############

# split fasta file where Ns occur

rule SPLIT_AT_Ns:
    input:
        "RESULTS/GENOME_ASSEMBLY/HIFIASM/{prefix}.p_ctgl2.fasta"
    output:
        "RESULTS/GENOME_ASSEMBLY/HIFIASM/{prefix}.p_ctgl2_split.fasta"
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/{prefix}.split_fa.log"
    params:
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"]
    shell: """
        split_fa {input} > {output}
    """

# xasm5=asm-to-ref mapping, for ~0.1 sequence divergence (assembly to self)
# -I=split index for every ~NUM input bases

rule MINIMAP2_GENOME_TO_SELF:
    input:
        "RESULTS/GENOME_ASSEMBLY/HIFIASM/{prefix}.p_ctgl2_split.fasta"
    output:
        "RESULTS/GENOME_ASSEMBLY/MINIMAP2/{prefix}.p_ctgl2.genome_split.paf"
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/{prefix}.minimap2a.log"
    params:
        chunk=config["CHUNKS"],
        i_split=config["INDEX_SPLIT"]
    shell: """
        minimap2 -I {params.i_split} -t {threads} -xasm5 -DP {input} {input} > {output}
    """

rule MINIMAP2_GENOME_TO_READS:
    input:
        splitf="RESULTS/GENOME_ASSEMBLY/HIFIASM/{prefix}.p_ctgl2_split.fasta",
        ccs_reads="DATA/{prefix}.fastq.gz"
    output:
        "RESULTS/GENOME_ASSEMBLY/MINIMAP2/{prefix}.p_ctgl2.reads_split.paf"
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/{prefix}.minimap2b.log"
    params:
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"],
        i_split=config["INDEX_SPLIT"]
    shell: """
        minimap2 -I {params.i_split} -t {threads} -x map-pb {input.splitf} {input.ccs_reads} > {output}
    """

#Calculate haploid/diploid coverage threshold and remove haplotype duplicates from assembly

rule PBCSTAT:
    input:
        "RESULTS/GENOME_ASSEMBLY/MINIMAP2/{prefix}.p_ctgl2.reads_split.paf"
    output:
        basecov="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}/PB.base.cov",
        stat="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}/PB.stat" #, PB.cov.wig
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/{prefix}.pbcstat.log"
    params:
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"],
        i_split=config["INDEX_SPLIT"],
        prefix="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}"
    shell: """
        pbcstat -O {params.prefix} {input}
    """

rule CALCUTS:
    input:
        "RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}/PB.stat"
    output:
        "RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}/cutoffs"
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/{prefix}.calcuts.log"
    params:
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"],
        i_split=config["INDEX_SPLIT"],
        prefix="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}"
    shell: """
        calcuts {input} > {output}
    """

rule PURGE_DUPS:
    input:
        basecov="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}/PB.base.cov",
        cutoffs="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}/cutoffs",
        genome_paf="RESULTS/GENOME_ASSEMBLY/MINIMAP2/{prefix}.p_ctgl2.genome_split.paf"
    output:
        "RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}_dups.bed"
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/{prefix}.purge_dups.log"
    params:
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"],
        i_split=config["INDEX_SPLIT"],
        prefix="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}"
    shell: """
        purge_dups -2 -c {input.basecov} -T {input.cutoffs} {input.genome_paf} > {output}
    """

rule GET_SEQS:
    input:
        purge_duped="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}_dups.bed",
        fasta="RESULTS/GENOME_ASSEMBLY/HIFIASM/{prefix}.p_ctgl2.fasta"
    output:
        "RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}_seqs_purged.hap.fa"
    threads:
        workflow.cores
    log:
        "RESULTS/GENOME_ASSEMBLY/LOG/{prefix}.get_seqs.log"
    params:
        loglevel=config["LOGLEVEL"],
        chunk=config["CHUNKS"],
        i_split=config["INDEX_SPLIT"],
        prefix="RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}_seqs_purged"
    shell: """
        get_seqs -e -p {params.prefix} {input.purge_duped} {input.fasta}
    """

rule INDEX_SEQS:
    input:
        "RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}_seqs_purged.hap.fa"
    output:
        "RESULTS/GENOME_ASSEMBLY/PURGE_DUPS/{prefix}_seqs_purged.hap.fa.fai"
    shell: """
        samtools faidx {input}
    """
