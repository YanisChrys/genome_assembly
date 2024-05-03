# add pretext here

#"RESULTS/HIC/BAM4VIZ/{basename}.bam"


#scaf_genom="RESULTS/HIC/YAHS/" + PREFIX + "_yahs_scaffolds_final.fa",
rule hicstuff_pipeline:
    input:
        scaf_genom=lambda wildcards: get_input(wildcards, SCAF_FILES, SCAF_BASENAMES, 'scaff_basename'),
        hic1=config["hic1"],
        hic2=config["hic2"]
    output:
        "RESULTS/HIC_CONTACT_MAP/{scaff_basename}/abs_fragments_contacts_weighted.cool",
        "RESULTS/HIC_CONTACT_MAP/{scaff_basename}/fragments_list.txt"
    threads:
        workflow.cores
    conda:
        "../envs/hicstuff.yaml"
    params:
        outdir="RESULTS/HIC_CONTACT_MAP/{scaff_basename}"
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

rule hicstuff_view_no_binning:
    input:
        frags="RESULTS/HIC_CONTACT_MAP/{scaff_basename}/fragments_list.txt",
        matrix="RESULTS/HIC_CONTACT_MAP/{scaff_basename}/abs_fragments_contacts_weighted.cool"
    output:
        "RESULTS/HIC_CONTACT_MAP/{scaff_basename}/{scaff_basename}_map.png"
    conda:
        "../envs/hicstuff.yaml"
    shell: """
        hicstuff view \
        --normalize \
        --frags ${input.frags} ${input.matrix} \
        --output {output}
    """

rule hicstuff_view_5k:
    input:
        frags="RESULTS/HIC_CONTACT_MAP/{scaff_basename}/fragments_list.txt",
        matrix="RESULTS/HIC_CONTACT_MAP/{scaff_basename}/abs_fragments_contacts_weighted.cool"
    output:
        "RESULTS/HIC_CONTACT_MAP/{scaff_basename}/{scaff_basename}_5kmap.png"
    conda:
        "../envs/hicstuff.yaml"
    shell: """
        hicstuff view \
        --binning 5kb \
        --normalize \
        --frags ${input.frags} ${input.matrix} \
        --output {output}
    """

rule hicstuff_view_10k:
    input:
        frags="RESULTS/HIC_CONTACT_MAP/{scaff_basename}/fragments_list.txt",
        matrix="RESULTS/HIC_CONTACT_MAP/{scaff_basename}/abs_fragments_contacts_weighted.cool"
    output:
        "RESULTS/HIC_CONTACT_MAP/{scaff_basename}/{scaff_basename}_10kmap.png"
    conda:
        "../envs/hicstuff.yaml"
    shell: """
        hicstuff view \
        --binning 10kb \
        --normalize \
        --frags ${input.frags} ${input.matrix} \
        --output {output}
    """

rule bandage:
    input:
        scaf_genom=lambda wildcards: get_input(wildcards, SCAF_FILES, SCAF_BASENAMES, 'scaff_basename')
    output:
        gfa_reduced=temp("RESULTS/BANDAGE/{scaff_basename}/{scaff_basename}_reduced.gfa"),
        gfa=temp("RESULTS/BANDAGE/{scaff_basename}/{scaff_basename}.gfa"),
        image="RESULTS/BANDAGE/{scaff_basename}/{scaff_basename}.graph.svg"
    threads:
        min(workflow.cores,3)
    conda:
        "../envs/bandage.yaml"
    params:
        outdir="RESULTS/BANDAGE/{scaff_basename}"
    conda:
        "../envs/bandage.yaml"
    shell: """
        gfastats {input} -o gfa > {output.gfa}

        # reduce size of graph
        # --scope depthrange
        # mindepth default: 10
        # maxdepth default: 100
        # distance is the number of node steps away to draw for the aroundnodes and aroundblast scopes: default 0
        Bandage reduce {output.gfa} {output.gfa_reduced} --mindepth 500.0 --maxdepth 1000000 --distance 3

        # produce assembly image
        # --colour uniform
        Bandage image {output.gfa_reduced} {output.image} --nodewidth 50 --depwidth 1
    """
