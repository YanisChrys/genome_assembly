############  GENOME ASSEMBLY WITH HIFIASM  ############

#-D FLOAT:    drop k-mers occurring >FLOAT*coverage times [5.0]
#-N INT  :     consider up to max(-D*coverage,-N) overlaps for each oriented read [100]

if INCLUDE_HIC:
    ruleorder:  hifiasm_hic > hifiasm_no_hic >  edit_hifiasm_output
elif not INCLUDE_HIC:
    ruleorder:  hifiasm_no_hic > hifiasm_hic > edit_hifiasm_output


rule hifiasm_hic:
    input:
        hifi_reads="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz",
        hic1=HIC1,
        hic2=HIC2
    output:
        "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.hic.p_ctg.gfa",
        "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.hic.a_ctg.gfa"
    threads:
        min(workflow.cores, 20)
    params:
        prefix="RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}",
        ops=lambda wc: hifiasm_options(wc.basename),
        hic=lambda wildcards, input: "--h1 " + input.hic1 + " --h2 " + input.hic2
        # hic=lambda input: "--h1 {} --h2 {}".format(input.hic1, input.hic2)
        # hic="--h1 {input.hic1} --h2 {input.hic2}"
    priority: 1
    conda:
        "../envs/hifiasm.yaml"
    shell: """
        hifiasm {params.ops} -t {threads} --primary {params.hic} -o {params.prefix} {input.hifi_reads}
    """

rule hifiasm_no_hic:
    input:
        hifi_reads="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"
    output:
        "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.p_ctg.gfa",
        "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.a_ctg.gfa"
    threads:
        min(workflow.cores, 20)
    params:
        prefix="RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}",
        ops=lambda wc: hifiasm_options(wc.basename),
        hic=""
    priority: 1
    conda:
        "../envs/hifiasm.yaml"
    shell: """
        hifiasm {params.ops} -t {threads} --primary {params.hic} -o {params.prefix} {input.hifi_reads}
    """

# create primary and alternate contigs for "-l2" by folding the file and writing as a fasta file
rule edit_hifiasm_output:
    input:
        hifigfa1=lambda wildcards: "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.hic.p_ctg.gfa" if INCLUDE_HIC else "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.p_ctg.gfa",
        hifigfa2=lambda wildcards: "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.hic.a_ctg.gfa" if INCLUDE_HIC else "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}.a_ctg.gfa"
    output:
        fasta1="RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}_primary.fasta",
        fasta2="RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_{basename}_alternate.fasta"
    threads:
        min(workflow.cores,2)
    priority: 1
    shell: """
        awk '/^S/{{print ">"$2;print $3}}' {input.hifigfa1} | fold > {output.fasta1}
        awk '/^S/{{print ">"$2;print $3}}' {input.hifigfa2} | fold > {output.fasta2}
    """

#gfastats testFiles/random2.gfa2.gfa -o fa // converts gfa to fasta