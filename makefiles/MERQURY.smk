
rule Merqury_hifiasm:
    input:
        meryldb="RESULTS/STATISTICS/FASTK/" + PREFIX + ".meryl",
        asm1=lambda wildcards: get_input_hifiasm(wildcards, HIFIASM_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[0],
        asm2=lambda wildcards: get_input_hifiasm(wildcards, HIFIASM_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[1]
    output:
        "RESULTS/STATISTICS/FASTK/{basename}_hifiasm.qv"
    threads:
        min(workflow.cores,15)
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        read_prefix=PREFIX,
        merqury_prefix="{basename}_hifiasm"
    shell: """
        cd RESULTS/STATISTICS/FASTK/
        export OMP_NUM_THREADS={threads}
        merqury.sh {input.meryldb} ../../../{input.asm1} ../../../{input.asm2} {params.merqury_prefix}
    """


rule Merqury_purge_dups:
    input:
        meryldb="RESULTS/STATISTICS/FASTK/" + PREFIX + ".meryl",
        asm1=lambda wildcards: get_input_purgedups(wildcards, PURGE_DUP_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[0],
        asm2=lambda wildcards: get_input_purgedups(wildcards, PURGE_DUP_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[1]
    output:
        "RESULTS/STATISTICS/FASTK/{basename}_purge_dups.qv"
    threads:
        min(workflow.cores,15)
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        read_prefix=PREFIX,
        merqury_prefix="{basename}_purge_dups"
    shell: """
        cd RESULTS/STATISTICS/FASTK/
        export OMP_NUM_THREADS={threads}
        merqury.sh {input.meryldb} ../../../{input.asm1} ../../../{input.asm2} {params.merqury_prefix}
    """

rule Merqury_scaffold:
    input:
        meryldb="RESULTS/STATISTICS/FASTK/" + PREFIX + ".meryl",
        asm1=lambda wildcards: get_input_scaff(wildcards, SCAF_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[0],
        asm2=lambda wildcards: get_input_scaff(wildcards, SCAF_FILES, HIFIASM_MERQURY_BASENAMES, 'basename')[1]
    output:
        "RESULTS/STATISTICS/FASTK/{basename}_yahs.qv"
    threads:
        min(workflow.cores,15)
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        read_prefix=PREFIX,
        merqury_prefix="{basename}_yahs"
    shell: """
        cd RESULTS/STATISTICS/FASTK/
        export OMP_NUM_THREADS={threads}
        merqury.sh {input.meryldb} ../../../{input.asm1} ../../../{input.asm2} {params.merqury_prefix}
    """
