#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q large.q,medium.q
#$ -N nextpolish
#$ -pe smp 51
#$ -e /share/pool/CompGenomVert/EpiAus/workflow/polishing/logfiles
#$ -o /share/pool/CompGenomVert/EpiAus/workflow/polishing/logfiles
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas

module load anaconda3/2022.05
conda activate genofish

export THREADS=$(expr ${NSLOTS} - 1)

snakemake \
    --snakefile makefiles/POLISH_NEXTPOLISH.smk \
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
    --rerun-incomplete \
    --rerun-triggers mtime \
    --stats "./stats_nextpolish.json" 

