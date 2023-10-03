#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas
#$ -N subsample_sequeltools


# activate module
# qsub -pe smp 64 -q medium.q subsample_sequeltools.sh

module load anaconda3/2022.05
module load samtools/1.10
conda activate seqtools

./SequelTools.sh -t S -u /share/pool/CompGenomVert/RawData/DresdenPhoxinus22062022/m54345U_220303_154000.subreads.bam -T lr -R 0.05 -f b -n $NSLOTS -v

conda deactivate
