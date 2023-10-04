#!/bin/bash

module load anaconda3/2022.05
conda activate genome_assembly

snakemake \
    --snakefile snakemake_sgehic.smk \
    --keep-going \
    --latency-wait 300 \
    -j ${THREADS} \
    --default-resources "tmpdir='/path/to/tmpdir'" \
    --verbose \
    --use-conda \
    --use-envmodules \
    --printshellcmds \
    --reason \
    --nolock \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --stats "./stats.json" \
    --report "./report.html"

