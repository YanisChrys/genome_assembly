#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q large.q,medium.q
#$ -N epiaus_flye
#$ -pe smp 21
#$ -e /share/pool/CompGenomVert/EpiAus/logfiles
#$ -o /share/pool/CompGenomVert/EpiAus/logfiles
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas

conda install -y -n flye -c conda-forge -c bioconda -c defaults "flye=2.9.2"
conda activate flye

set -xe 
export THREADS=$(expr ${NSLOTS} - 1)

# Ensure necessary directories and files exist
outdir=${PWD}/"RESULTS/GENOME_ASSEMBLY_FLYE"
reads="/share/pool/CompGenomVert/EpiAus/workflow/DATA/DECONTAMINATED/EpiAus_202401_3SMRT_kraken_unclassified.fastq.gz"
genome_size="1869647118"

mkdir -p ${outdir}

flye --threads $THREADS \
    --out-dir $outdir \
    --genome-size $genome_size \
    --pacbio-hifi $reads
    #--scaffold
    #--keep-haplotypes
