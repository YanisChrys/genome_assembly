############  GENOME ASSEMBLY WITH HIFIASM  ############

# HiFiasm
# also try -l0 and -l3

# l: 0-3 level of purging
# primary: output a primary assembly and an alternate assembly
# o: prefix of output files
# t: threads

rule HIFIASM:
    input:
        "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq"
    output:
        "RESULTS/GENOME_ASSEMBLY/" + PREFIX + ".p_ctg.gfa",
        "RESULTS/GENOME_ASSEMBLY/" + PREFIX + ".a_ctg.gfa"
    threads:
        workflow.cores
    params:
        prefix="RESULTS/GENOME_ASSEMBLY/" + PREFIX,
        chunk=config["CHUNKS"]
    conda:
        "../envs/hifiasm.yaml"
    shell: """
        hifiasm -l2 -t {threads} --primary -o {params.prefix} {input}
    """

# create primary and alternate contigs for "-l2" by folding the file and writing as a fasta file
rule EDIT_HIFIASM_OUTPUT:
    input:
        hifigfa1="RESULTS/GENOME_ASSEMBLY/" + PREFIX + ".p_ctg.gfa",
        hifigfa2="RESULTS/GENOME_ASSEMBLY/" + PREFIX + ".a_ctg.gfa"
    output:
        fasta1="RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_primary.fasta",
        fasta2="RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_alternate.fasta"
    shell: """
        awk '/^S/{{print ">"$2;print $3}}' {input.hifigfa1} | fold > {output.fasta1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hifigfa2} | fold > {output.fasta2}
    """

#gfastats testFiles/random2.gfa2.gfa -o fa // converts gfa to fasta