


# next we pay special attention to some functions that are costly to run and should therefore
# only run on the terminal assembly, if at all:
TERMINAL_ASSEMBLY = []
if DO_DEDUP and not DO_HIC:
    TERMINAL_ASSEMBLY = expand("RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa", basename=HIFIASM_BASENAMES)
if DO_DEDUP and DO_HIC:
    TERMINAL_ASSEMBLY = SCAF_FILES
TERMINAL_BASENAMES = [os.path.basename(f).split('.')[0] for f in TERMINAL_ASSEMBLY]

    rule all: 
        input:
            # read mapping
            qualimap_hifi = expand("RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/{basename}_report.pdf", basename = TERMINAL_BASENAMES)
            qualimap_hic = expand("RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIC/{basename}_report.pdf", basename = SCAF_BASENAMES)
            fastq = expand("RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/{basename}.fastq", basename = TERMINAL_BASENAMES)
            flagstat = expand("RESULTS/STATISTICS/READ_MAPPING/{basename}.flagstat", basename = TERMINAL_BASENAMES)
            scafmap = expand("RESULTS/HIC/BAM4VIZ/{basename}.bam", basename = SCAF_BASENAMES)
            qualimap_hifi + qualimap_hic + fastq + flagstat + scafmap


# qualimap

rule map_hifi_to_db:
    input:
        ref=lambda wildcards: get_input(wildcards, TERMINAL_ASSEMBLY, TERMINAL_BASENAMES, 'basename'),
        reads = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz",
    output:
        "RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/{basename}.bam"
    threads:
        workflow.cores
    params:
        rg_group="'@RG\\tSM:{basename}\\tID:{basename}\\tPL:PACBIO'"
    conda:
        "../envs/samtools_minimap.yaml"
    shell: """
        minimap2  -a --MD --eqx -x map-hifi  \
        -k 20 \
        -R {params.rg_group} \
        -t {threads} \
        {input.ref} {input.reads} | samtools view -Sbh  | samtools sort -o {output} -
    """

rule map_scaf2hic:
    input:
        fasta=lambda wildcards: get_input(wildcards, TERMINAL_ASSEMBLY, TERMINAL_BASENAMES, 'basename'),
        hic1=HIC1,
        hic2=HIC2
    output:
        temp("RESULTS/HIC/BAM4VIZ/{basename}.bam")
    threads:
        workflow.cores
    conda:
        "../envs/bwa.yaml"
    params:
        rggroup = "'@RG\tID:{basename}\tSM:{basename}\tLB:{basename}\tPL:ILLUMINA'"
    shell: """
        bwa-mem2 index {input.fasta}
        bwa-mem2 mem -t{threads} -B8 -M -R {params.rggroup} {input.fasta} {input.hic1} {input.hic2} | samtools view -@ {threads} -bS -o - > {output}
    """


# seqkit seq -Q64 -V  Read1.fastq > Read1.sanger.fastq

# 5) extract good mapped reads
# exclude:
# secondary alignments (256)
# supplementary alignments (2048)
# unmapped (4)
# 3844 = (read unmapped, 
    # read fails platform/vendor quality checks,
    # not primary alignment,
    # PCR or optical duplicate,
    # supplementary alignment)

rule extract_only_mapped_reads:
    input:
        bam="RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/{basename}.bam"
    output:
        fastq=temp("RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/{basename}.fastq"),
        flagstat="RESULTS/STATISTICS/READ_MAPPING/{basename}.flagstat"
    conda:
        "../envs/samtools_minimap.yaml"
    threads:
        workflow.cores
    shell: """
        samtools fastq -@ {threads} -F 256 -F 2048 -F 4 {input.bam} > {output.fastq}
        samtools flagstat -O tsv {output.fastq} > {output.flagstat}
    """

rule qualimap_hifi:
    input:
        "RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/{basename}.bam"
    output:
        qualimap="RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/{basename}_report.pdf"
    conda:
        "../envs/qualimap.yaml"
    threads:
        min(workflow.cores,20)
    params:
        dir="RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIFI/",
        report="{basename}_report.pdf"
    shell: """
        qualimap bamqc -bam {input} -c -outdir {params.dir} -outfile {params.report} -nt {threads} --java-mem-size=10G
    """

rule qualimap_hic:
    input:
        "RESULTS/HIC/BAM4VIZ/{basename}.bam"
    output:
        qualimap="RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIC/{basename}_report.pdf"
    conda:
        "../envs/qualimap.yaml"
    threads:
        min(workflow.cores,20)
    params:
        dir="RESULTS/STATISTICS/READ_MAPPING/QUALIMAP_HIC/",
        report="{basename}_report.pdf"
    shell: """
        qualimap bamqc -bam {input} -c -outdir {params.dir} -outfile {params.report} -nt {threads} --java-mem-size=10G
    """



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
