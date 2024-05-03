
#possible options to try:
# --scaffold
#--keep-haplotypes

rule format_genome_size_for_flye:
    input:
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_genome_size.txt"
    output:
        formatted_genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_formatted_genome_size.txt"
    run:
        def format_genome_size(size):
            size = int(size)
            if size >= 1000000000:
                return f'{size / 1000000000:.1f}g'
            elif size >= 1000000:
                return f'{size / 1000000:.1f}m'
            elif size >= 1000:
                return f'{size / 1000:.1f}k'
            else:
                return f'{size}'

        with open(input.genome_size_file, 'r') as file:
            raw_size = file.read().strip()
        formatted_size = format_genome_size(raw_size)
        with open(output.formatted_genome_size_file, 'w') as file:
            file.write(formatted_size)

rule flye:
    input:
        reads = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz",
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_formatted_genome_size.txt"
    output:
        outdir = directory("RESULTS/GENOME_ASSEMBLY_FLYE/"),
        asm = "RESULTS/GENOME_ASSEMBLY_FLYE/assembly.fasta"
    log:
        "RESULTS/LOGS/FLYE.log"
    threads:
        min(workflow.cores,20)
    conda:
        "../envs/flye.yaml"
    params:
        genome_size=lambda wildcards, input: open(input.genome_size_file).read().strip()
    shell: """
        (flye \
        --threads {threads} \
        --out-dir {output.outdir} \
        --genome-size {params.genome_size} \
        --pacbio-hifi {input.reads}) 2> {log}
    """