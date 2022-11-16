############  GENOME ASSEMBLY WITH HIFIASM  ############

# HiFiasm
# also try -l0 and -l3

# l: 0-3 level of purging
# primary: output a primary assembly and an alternate assembly
# o: prefix of output files
# t: threads

# inspired by:
# https://snakemake-wrappers.readthedocs.io/en/latest/wrappers/hifiasm.html
# Input: list of files
# output: list of files that differ by their extension
# run hifiasm so that if hic is provided it is used (else [" {hic1} {hic2}"] string will be empty)
rule RUN_HIFIASM:
    input:
        "DATA/{prefix}.fastq.gz"
    output:
        "RESULTS/GENOME_ASSEMBLY/{prefix}.p_ctg.gfa",
        "RESULTS/GENOME_ASSEMBLY/{prefix}.a_ctg.gfa"
        # multiext("RESULTS/GENOME_ASSEMBLY/{prefix}.",
        #     "a_ctg.gfa",
        #     "a_ctg.lowQ.bed",
        #     "a_ctg.noseq.gfa",
        #     "p_ctg.gfa",
        #     "p_ctg.lowQ.bed",
        #     "p_ctg.noseq.gfa",
        #     "p_utg.gfa",
        #     "p_utg.lowQ.bed",
        #     "p_utg.noseq.gfa",
        #     "r_utg.gfa",
        #     "r_utg.lowQ.bed",
        #     "r_utg.noseq.gfa",
        # )
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{prefix}.hifiasm.log"
    params:
        prefix="RESULTS/GENOME_ASSEMBLY/{prefix}",
        chunk=config["CHUNKS"]
    envmodules:
        "hifiasm/0.16.1"
    shell: """
        hifiasm -l2 -t {threads} --primary -o {params.prefix} {input}
    """

# create primary and alternate contigs for "-l2" by folding the file and writing as a fasta file
rule EDIT_HIFIASM_OUTPUT:
    input:
        hifigfa1="RESULTS/GENOME_ASSEMBLY/{prefix}.p_ctg.gfa",
        hifigfa2="RESULTS/GENOME_ASSEMBLY/{prefix}.a_ctg.gfa"
    output:
        fasta1="RESULTS/GENOME_ASSEMBLY/{prefix}_primary.fasta",
        fasta2="RESULTS/GENOME_ASSEMBLY/{prefix}_alternate.fasta"
    shell: """
        awk '/^S/{{print ">"$2;print $3}}' {input.hifigfa1} | fold > {output.fasta1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hifigfa2} | fold > {output.fasta2}
    """
