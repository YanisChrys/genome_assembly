



rule run_hypo:
    input:
        asm="RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fa"
    output:
        out1=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_hypo-1.fa"),
        out2=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_hypo-2.fa"),
        out3="RESULTS/ASM_POLISHING_HYPO/{sample}_hypo-3.fa",
        bam1=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_hypo-1.bam"),
        bam2=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_hypo-2.bam"),
        bam3=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_hypo-3.bam")
    threads: 
        workflow.cores
    params:
        genome_size=config["genome_size"],
        coverage = lambda wildcards: calculate_cov(config["total_read_length"],config["genome_size"]),
        pacbioreads=config["reads"],
        fofn="RESULTS/ASM_POLISHING_HYPO/inputs.fofn"
    conda:
        "../envs/hypo.yaml"
    shell: """
        cat <<EOF > {params.fofn} 
{params.pacbioreads}
EOF

        minimap2 --secondary=no --MD --eqx -ax map-hifi -t {threads} {input.asm} $(echo {input.pacbioreads}) | samtools view -Sb - | samtools sort -o {output.bam1} -@ {threads} -
        samtools index {output.bam1}
        hypo -d {input.asm} -r @{output.fofn} -s {params.genome_size} -c {params.coverage} -b {output.bam1} -p 80 -t {threads} -o {output.out1}

        minimap2 --secondary=no --MD --eqx -ax map-hifi -t {threads} {output.out1} $(echo {input.pacbioreads}) | samtools view -Sb - | samtools sort -o {output.bam2} -@ {threads} -
        samtools index {output.bam2}
        hypo -d {output.out1} -r @{output.fofn} -s {params.genome_size} -c {params.coverage} -b {output.bam2} -p 80 -t {threads} -o {output.out2}

        minimap2 --secondary=no --MD --eqx -ax map-hifi -t {threads} {output.out2} $(echo {input.pacbioreads}) | samtools view -Sb - | samtools sort -o {output.bam3} -@ {threads} -
        samtools index {output.bam3}
        hypo -d {output.out2} -r @{output.fofn} -s {params.genome_size} -c {params.coverage} -b {output.bam3} -p 80 -t {threads} -o {output.out3}
"""

rule align_reads:
    input:
        ref="{{path}}/{sample}_{{treatment}}.fa"
    output:
        bam="{{path}}/ASM_POLISHING_DEEPVARIANT/alignments/{sample}_{{treatment}}.bam"
    threads:
        min(workflow.cores,20)
    conda:
        "../envs/pbmm2.yaml"
    params:
        pbmm2_ops=config["pbmm2"]["ops"],
        reads=config["reads"]
    shell: """
        pbmm2 align {input.ref} \
        {params.reads} {output.bam} \
        {params.pbmm2_ops} \
        -j {threads} \
        --rg '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}'
    """



# split genome into chromosomes
# two-pass: read file twice to reduce memory usage
# split by the full id name of each chromosome
checkpoint split_by_chromosomes:
    input:
        "{{path}}/{sample}_{{treatment}}.fa"
    output:
        dynamic("{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.fasta")
    conda:
        "../envs/seqkit.yaml"
    params:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/"
    shell: """
        seqkit split --quiet --two-pass -i --by-id-prefix "" \
        --out-dir {params} {input}
        """

rule index_chromosomes:
    input:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.fasta"
    output:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.fasta.fai"
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools faidx {input}
    """

rule create_bed:
    input:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.fasta.fai"
    output:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.bed"
    shell: """
        awk '{{printf("%s\\t0\\t%s\\n",$1,$2)}}' {input}  > {output}
    """

rule split_aligned_reads:
    input:
        bam="{{path}}/ASM_POLISHING_DEEPVARIANT/alignments/{sample}_{{treatment}}.bam",
        bed="{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.bed"
    output:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/alignments_split/{sample}_{{treatment}}/{chrom}_filt.bam"
    threads:
        min(workflow.cores,20)
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools view -@{threads} -bh -F 2308 -M -L {input.bed} {input.bam} > {output}
        samtools index {output}
    """

rule run_deepvariant:
    input:
        ref = "{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.fasta",
        bam = "{{path}}/ASM_POLISHING_DEEPVARIANT/alignments_split/{sample}_{{treatment}}/{chrom}_filt.bam"
    output:
        vcf = "{{path}}/ASM_POLISHING_DEEPVARIANT/deepvariant/{sample}_{{treatment}}/{chrom}/all.filt.dv.vcf"
    threads: 
        min(workflow.cores,20)
    container:
        "/share/scientific_bin/singularity/containers/deepvariant_1.2.0.sif"
    params:
        "output/intermediate_files/"
    shell: """ 
        run_deepvariant \
        --model_type=PACBIO \
        --ref={input.ref} \
        --reads={input.bam} \
        --output_vcf={output.vcf} \
        --intermediate_results_dir {params} \
        --num_shards={threads}
    """

rule filter_deepvariant:
    input:
        chrom="{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.fasta",
        vcf="{{path}}/ASM_POLISHING_DEEPVARIANT/deepvariant/{sample}_{{treatment}}/{chrom}/all.filt.dv.vcf"
    output:
        vcf="{{path}}/ASM_POLISHING_DEEPVARIANT/deepvariant_filtered/{sample}_{{treatment}}/{chrom}/all.filt.dv.filt.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell: """
        bcftools view -f 'PASS' -i 'GT="1/1"' {input.vcf} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
    """

rule filter_per_chromosome:
    input:
        fasta="{{path}}/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_{{treatment}}/{chrom}.fasta",
        vcf="{{path}}/ASM_POLISHING_DEEPVARIANT/deepvariant_filtered/{sample}_{{treatment}}/{chrom}/all.filt.dv.filt.vcf.gz"
    output:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/chrom_consensus/{sample}_{{treatment}}/{chrom}.consensus.fasta"
    conda:
        "../envs/bcftools.yaml"
    shell: """
        bcftools consensus -f {input.fasta} {input.vcf} > {output}
    """



rule merge_chroms:
    input:
        dynamic("{{path}}/ASM_POLISHING_DEEPVARIANT/chrom_consensus/{sample}_{{treatment}}/{chrom}.consensus.fasta")
    output:
        "{{path}}/ASM_POLISHING_DEEPVARIANT/{sample}_{{treatment}}_deepvariant.fasta"
    shell:
        "cat {input} > {output}"
