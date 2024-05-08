#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q large.q,medium.q
#$ -N stats_pri
#$ -pe smp 31
#$ -e /share/pool/CompGenomVert/EpiAus/workflow/logfiles
#$ -o /share/pool/CompGenomVert/EpiAus/workflow/logfiles
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas

module load anaconda3/2022.05
conda activate genofish

export THREADS=$(expr ${NSLOTS} - 1)

snakemake \
    --snakefile STATS.smk \
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
    --rerun-incomplete --until MerquryFK 
    
    
    
    #\
    #--stats "./stats.json" 


# snakemake \
#     --snakefile VISUALIZE.smk \
#     --keep-going \
#     --latency-wait 300 \
#     -j ${THREADS} \
#     --default-resources "tmpdir='/path/to/tmpdir'" \
#     --verbose \
#     --use-conda \
#     --use-envmodules \
#     --printshellcmds \
#     --reason \
#     --nolock \
#     --rerun-incomplete -n \
#     --stats "./stats.json" 

    #--report "./report.html"

