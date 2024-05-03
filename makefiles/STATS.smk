############  INITIALIZE PRIMARY STATISTICAL MODULE AND TABULAR ARRAY ############

# Pipeline Programs and purpose:
# - SEQKIT: Used for initial read statistics from the CCS reads in rules SEQKIT_CCS_STATS, SEQKIT_ASSEMBLY_STATS, and SEQKIT_YAHS_STATS.
# - FastK: Used for k-mer analysis in rule FASTK_TABLE.
# - Histex: Used for simplifying FastK hist table in rule HISTEX.
# - GeneScopeFK: Used for generating plots describing genome properties in rule GENESCOPE_FK.
# - KatGC: Used for generating 3D heat map or contour map of the frequency of a k-mer versus its' GC content in rule KATGC.
# - MerquryFK: Used for evaluating the completeness of the assembly in rule MerquryFK.
# - PloidyPlot: Used for generating ploidy plots in rule PLOIDYPLOT.
# - QUAST: Used for evaluating assembly quality in rule quast.
# - GT SEQSTATS: Used for generating sequence statistics in rule gtseqstat.
# - GFASTATS: Used for generating GFA statistics in rule gfastats.

rule seqkit_hifi_stats:
    input:
        "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"
    output:
        "RESULTS/STATISTICS/SEQKIT/" + PREFIX + ".ccs.readStats.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """

rule fastq_hic_stats:
    input:
        "RESULTS/HIC/YAHS/" + PREFIX + "_yahs_scaffolds_final.fa"
    output:
        "RESULTS/STATISTICS/SEQKIT/" + PREFIX + "._yahs_scaffolds_final_assemblyStat.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
"""

# FastK
# kmers are k-sized string chunks - here it is DNA
# FastK splits the DNA into kmers and counts how many timers it sees each one
# need to provide absolute path of the module's locations
# See the README.mk for the prerequisite steps to running the programs from Fastk and Merqury packages
rule meryl:
    input:
        "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"
    output:
        directory("RESULTS/STATISTICS/FASTK/" + PREFIX + ".meryl")
    threads:
        min(workflow.cores,15)
    priority: 1
    params:
        wd="RESULTS/STATISTICS/FASTK/",
        out=PREFIX + ".meryl",
        kmers=config["modules"]["meryl"]["kmer"]
    conda:
        "../envs/meryl.yaml"
    shell: """
        cd {params.wd}
        meryl count k={params.kmers} ../../../{input} threads={threads} memory=200g output {params.out}
    """

############  GENOME SIZE ESTIMATE WITH GENESCOPE.FK  ############
rule meryl_hist:
    input:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + ".meryl"
    output:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist"
    conda:
        "../envs/meryl.yaml"
    threads:
        1
    shell:
        "meryl histogram {input} > {output}"

rule genomescope:
    input:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist"
    output:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + "_genomescope/summary.txt"
    conda:
        "../envs/genomescope.yaml"
    threads:
        1
    params:
        outdir = "RESULTS/STATISTICS/FASTK/" + PREFIX + "_genomescope",
        ploidy=2,
        kmer=config["modules"]["meryl"]["kmer"]
    shell:
        "genomescope2 -p {params.ploidy} -k {params.kmer} -i {input} -o {params.outdir}"



############  READ ANALYSIS  ############


rule smudgeplot:
    input:
        db="RESULTS/STATISTICS/FASTK/" + PREFIX + ".meryl",
        hist="RESULTS/STATISTICS/FASTK/" + PREFIX + ".hist"
    output:
        "RESULTS/STATISTICS/FASTK/" + PREFIX + ".png"
    conda:
        "../envs/smudgeplot.yaml"
    threads:
        min(workflow.cores,10)
    params:
        prefix=PREFIX,
        outdir="RESULTS/STATISTICS/FASTK"
    shell: """
        L=$(smudgeplot.py cutoff {input.hist} L)
        U=$(smudgeplot.py cutoff {input.hist} U)

        meryl print less-than 1500 greater-than 10 threads={threads} memory=10G {input.db} | sort | \
        smudgeplot.py hetkmers -o {params.outdir}/{params.prefix}_L${{L}}_U${{U}} 
        smudgeplot_plot.R -L ${{L}} -i {params.outdir}/{params.prefix}_L${{L}}_U${{U}}_coverages.tsv -o {params.outdir}/{params.prefix}
    """


# ASSEMBLY STATS

rule seqkit_assembly_stats:
    input:
        lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS/SEQKIT/{basename}.assemblyStat.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """

# QUAST

rule quast:
    input:
        contigs = lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS/QUAST/{basename}/quast.log"
    threads: 
        min(workflow.cores,20)
    priority: 1
    params:
        "RESULTS/STATISTICS/QUAST/{basename}"
    conda:
        "../envs/quast.yaml"
    shell: """
        quast.py \
        --eukaryote \
        --large \
        -t {threads} \
        --min-contig 300 \
        -o {params} \
        {input.contigs}
    """

# GT SEQSTATS

# extract genome size to a file to help snakemake know it needs to be calculated before running the stats
rule extractgenomeSize:
    input:
        genescope_output="RESULTS/STATISTICS/FASTK/" + PREFIX + "_genomescope/summary.txt"
    output:
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_genome_size.txt"
    params:
        genome_size = lambda wildcards, input: extract_max_genome_haploid_length(input.genescope_output)
    shell:
        "echo {params.genome_size} > {output.genome_size_file}"

rule gtseqstat:
    input:
        assembly= lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename'),
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_genome_size.txt"
    output:
        stats="RESULTS/STATISTICS/GTSEQSTATS/{basename}_gtseqstat.txt"
    priority: 1
    params:
        genome_size=lambda wildcards, input: open(input.genome_size_file).read().strip()
    conda:
        "../envs/genometools.yaml"
    shell: """
        # module load genometools/1.6.2
        gt seqstat -contigs -genome {params.genome_size} \
        {input.assembly} > {output.stats}
    """

# GFASTATS

rule gfastats:
    input:
        fasta=lambda wildcards: get_input(wildcards, ASSEMBLY_FILES, ASSEMBLY_BASENAMES, 'basename'),
        genome_size_file="RESULTS/STATISTICS/" + PREFIX + "_genome_size.txt"
    output:
        "RESULTS/STATISTICS/GFASTATS/{basename}_gfastats_report.txt"
    conda:
        "../envs/gfastas.yaml"
    threads:
        min(workflow.cores,10)
    priority: 1
    params:
        genome_size=lambda wildcards, input: str(int(open(input.genome_size_file).read().strip()))
    shell: """
        gfastats \
        {input.fasta} \
        {params.genome_size} \
        --threads {threads} \
        --tabular \
        --nstar-report \
        --seq-report  > {output}
    """

# --discover-paths \
# --out-size scaffolds \

rule qualimap:
    input:
        "RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG_deduped.bam"
    output:
        report="RESULTS/STATISTICS/QUALIMAP/{basename}/report.pdf"
    conda:
        "../envs/qualimap.yaml"
    threads:
        min(workflow.cores,10)
    params:
        dir="RESULTS/STATISTICS/QUALIMAP/{basename}/"
    shell: """
        qualimap bamqc -bam {input} -c -outdir {params.dir} -outformat pdf -nt {threads} --java-mem-size=10G
    """

# visualize

rule pretext:
    input:
        bam="RESULTS/HIC/COMBINED/{basename}_mapped_hic_reads_sorted_withRG_deduped.bam"
    output:
        pretext="RESULTS/STATISTICS/PRETEXT/{basename}.pretext"
    conda:
        "../envs/pretext.yaml"
    threads:
        min(workflow.cores,5)
    shell: """
        samtools view -h {input.bam} | PretextMap -o {output.pretext}
        PretextSnapshot -m {output.pretext} -f png -r 1000 -c 5 --sequences '=full' --minTexels 64 --gridSize 1 --gridColour black '' -o output --prefix contigs
    """

rule hicstuff_pipeline:
    input:
        scaf_genom=lambda wildcards: get_input(wildcards, SCAF_FILES, SCAF_BASENAMES, 'scaff_basename'),
        hic1=HIC1,
        hic2=HIC2
    output:
        "RESULTS/STATISTICS/HICSTUFF/{scaff_basename}/abs_fragments_contacts_weighted.cool",
        "RESULTS/STATISTICS/HICSTUFF/{scaff_basename}/fragments_list.txt"
    threads:
        min(workflow.cores,20)
    conda:
        "../envs/hicstuff.yaml"
    params:
        outdir="RESULTS/HICSTUFF/{scaff_basename}"
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
        frags="RESULTS/STATISTICS/HICSTUFF/{scaff_basename}/fragments_list.txt",
        matrix="RESULTS/STATISTICS/HICSTUFF/{scaff_basename}/abs_fragments_contacts_weighted.cool"
    output:
        "RESULTS/STATISTICS/HICSTUFF/{scaff_basename}/{scaff_basename}_map.png"        
    conda:
        "../envs/hicstuff.yaml"
    threads:
        min(workflow.cores,20)
    shell: """
        hicstuff view \
        --normalize \
        --frags ${input.frags} ${input.matrix} \
        --output {output}
    """
