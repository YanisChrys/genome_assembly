#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q large.q,medium.q
#$ -N pol_epid1
#$ -pe smp 41
#$ -e /share/pool/CompGenomVert/RefGenomes/EpiAus/workflow/polishing/logfiles
#$ -o /share/pool/CompGenomVert/RefGenomes/EpiAus/workflow/polishing/logfiles
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas

# h1 = hypo
# d1 = deepvariant
# d2 = deepvariant on hypo

module load anaconda3/2022.05
module load singularity/3.10.5_fix
conda activate snakemake8

mkdir -p $PWD/temp
export TMPDIR="$PWD/temp"
export THREADS=$(expr ${NSLOTS} - 1)


# snakemake \
#     --snakefile makefiles/POLISH_DEEPVARIANT_old_prep.smk \
#     --keep-going \
#     --latency-wait 300 \
#     -j ${THREADS} \
#     --default-resources "tmpdir='/share/pool/ychrysostomakis/tmp'" \
#     --verbose \
#     --use-conda \
#     --use-envmodules \
#     --printshellcmds \
#     --nolock \
#     --rerun-incomplete \
#     --rerun-triggers mtime

snakemake \
    --snakefile makefiles/POLISH_DEEPVARIANT_old.smk \
    --keep-going \
    --latency-wait 300 \
    -j ${THREADS} \
    --verbose \
    --use-conda \
    --use-envmodules \
    --use-singularity \
    --singularity-args "--home $PWD" \
    --singularity-args "--bind $TMPDIR:$TMPDIR" \
    --default-resources "tmpdir='$TMPDIR'" \
    --printshellcmds \
    --nolock \
    --rerun-incomplete \
    --rerun-triggers mtime


