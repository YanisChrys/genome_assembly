#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q medium.q
#$ -N genome_assembly
#$ -pe smp 60
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas

module load anaconda3/2022.05
conda activate genome_assembly

#one core will be used by snakemake to monitor the other processes
THREADS=$(expr ${NSLOTS} - 1)

snakemake \
    --snakefile snakemake_hic \
    --keep-going \
    --latency-wait 60 \
    --use-envmodules \
    --use-conda \
    --cores ${THREADS} \
    --verbose \
    --printshellcmds \
    --reason \
    --resources mem_mb=1000 \
    --nolock 
    
    #--rerun-triggers mtime

