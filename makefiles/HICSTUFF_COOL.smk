
#scaf_genom="RESULTS/HIC/YAHS/" + PREFIX + "_yahs_scaffolds_final.fa",
rule HICSTUFF_PIPELINE:
    input:
        scaf_genom=lambda wildcards: get_input(wildcards, SCAF_FILES, SCAF_BASENAMES, 'scaff_basename'),
        hic1=config["HIC1"],
        hic2=config["HIC2"]
    output:
        "RESULTS/HIC_CONTACT_MAP_{scaff_basename}/abs_fragments_contacts_weighted.cool",
        "RESULTS/HIC_CONTACT_MAP_{scaff_basename}/fragments_list.txt"
    threads:
        workflow.cores
    conda:
        "../envs/hicstuff.yaml"
    params:
        outdir="RESULTS/HIC_CONTACT_MAP_{scaff_basename}"
    shell: """
        hicstuff pipeline \
            -t {threads} \
            --matfmt="cool" \
            --plot \
            --duplicates \
            --distance-law \
            --genome {input.scaf_genom} \
            --outdir {params.outdir} \
            {input.hic1} \
            {input.hic2}
    """

rule HICSTUFF_VIEW_NO_BINNING:
    input:
        frags="RESULTS/HIC_CONTACT_MAP_{scaff_basename}/fragments_list.txt",
        matrix="RESULTS/HIC_CONTACT_MAP_{scaff_basename}/abs_fragments_contacts_weighted.cool"
    output:
        "RESULTS/HIC_CONTACT_MAP_{scaff_basename}/{scaff_basename}_map.png"
    conda:
        "../envs/hicstuff.yaml"
    shell: """
        hicstuff view \
        --normalize \
        --frags ${input.frags} ${input.matrix} \
        --output {output}
    """

rule HICSTUFF_VIEW_5k:
    input:
        frags="RESULTS/HIC_CONTACT_MAP_{scaff_basename}/fragments_list.txt",
        matrix="RESULTS/HIC_CONTACT_MAP_{scaff_basename}/abs_fragments_contacts_weighted.cool"
    output:
        "RESULTS/HIC_CONTACT_MAP_{scaff_basename}/{scaff_basename}_5kmap.png"
    conda:
        "../envs/hicstuff.yaml"
    shell: """
        hicstuff view \
        --binning 5kb \
        --normalize \
        --frags ${input.frags} ${input.matrix} \
        --output {output}
    """

rule HICSTUFF_VIEW_10k:
    input:
        frags="RESULTS/HIC_CONTACT_MAP_{scaff_basename}/fragments_list.txt",
        matrix="RESULTS/HIC_CONTACT_MAP_{scaff_basename}/abs_fragments_contacts_weighted.cool"
    output:
        "RESULTS/HIC_CONTACT_MAP_{scaff_basename}/{scaff_basename}_10kmap.png"
    conda:
        "../envs/hicstuff.yaml"
    shell: """
        hicstuff view \
        --binning 10kb \
        --normalize \
        --frags ${input.frags} ${input.matrix} \
        --output {output}
    """

