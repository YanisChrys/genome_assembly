#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q large.q
#$ -pe smp 100
#$ -N genome_assembly
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas

module load anaconda3/2022.05
conda activate genofish

#one core will be used by snakemake to monitor the other processes

THREADS=$(expr ${NSLOTS} - 1)

snakemake \
    --snakefile snakemake_sgehic.smk \
    --keep-going \
    --latency-wait 300 \
    -j ${THREADS} \
    --default-resources "tmpdir='/share/pool/ychrysostomakis/tmp'" \
    --verbose \
    --use-conda \
    --use-envmodules \
    --printshellcmds \
    --reason \
    --nolock \
    --rerun-triggers mtime \
    --rerun-incomplete --conda-create-envs-only --until HIFIASM  \
    --stats "./stats.json" \
    --report "./report.html"

