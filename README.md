# genome_assembly

## NOTE: Maybe run pipeline?

```
# salsa2.3 needs python 2.7 The syntax of Python 3 is different
conda create --name genome_assembly
conda activate genome_assembly
conda config --add channels conda-forge
conda config --add channels bioconda  
conda install -y purge_dups pbmm2 pbccs pbgcpp bedtools r-essentials r-argparse r-ggplot2 picard r-scales r-viridis samtools snakemake bcftools seqkit freebayes python=2.7 actc hifiasm bam2fastx bwa pairtools pairix r-base r-minpack.lm busco merqury openjdk=11
```

The following packages need to be set up separately by following the instructions on the github page or by installing them in whatever way is appropriate for your system:

[FastK](https://github.com/thegenemyers/FASTK) 

[MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK) 

[GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK) 

[Salsa2](https://github.com/VGP/vgp-assembly/tree/master/pipeline/salsa)

[DeepVariant](https://github.com/google/deepvariant)  

[higlass](https://github.com/higlass/higlass)  

### Activate conda environment with:
``` 
conda activate genome_assembly
```



## Install rule-specific conda environments
Inside the folder `/envs` are `.yaml` files used to install conda environments needed for some of the operations in the workflow. When planning to run the workflow offline, these environments need to be installed beforehand _by snakemake_ when you are online. To do that, you need to use the following command on the terminal (in a cluster system that is the head node of the user) :

```
 snakemake -s snakemake_hic --cores 5 -p -r -w 5 --verbose --conda-create-envs-only --use-conda --conda-frontend conda
```
This needs to be done every time the environment files are changed.

### Do a dry run on 5 cores with:
```
snakemake -s snakemake_hic --dry-run --cores 5 -p -r -w 5 --verbose
```


### Run on SGE cluster with 60 cores with:
```
qsub -pe smp 60 -q medium.q SGE_genome_assembly.sh
```

### Create graph of jobs with:
```
snakemake -s PATH/TO/SMKFILE/snakemake_hic --dag --forceall | dot -Tpdf > graph_of_jobs.pdf
```
## BUSCO:

For BUSCO, it is preferable to use it `offline`. To do that, download the desired dataset from [here](https://busco-data.ezlab.org/v5/data/lineages/), unpack and place in folder with:
```
curl -O <link>
tar -xf <file> # becomes <folder>
mkdir busco_downloads
mkdir busco_downloads/lineages
mv <folder> busco_downloads/lineages
```


# Samples:
/share/pool/CompGenomVert/RawData/DresdenPhoxinus22062022

## HiC:
fPhoPho_R001_S001_R1.fastq.gz  
fPhoPho_R001_S001_R2.fastq.gz  

### subreads:  
m54345U_220213_143410.sts.xml  
m54345U_220213_143410.subreads.bam  
m54345U_220213_143410.subreads.bam.pbi  
m54345U_220213_143410.subreadset.xml  

### enzyme motif:  
GATC,GANTC,CTNAG,TTAA  

### HiFi:  
m54345U_220213_143410.deepconsensusFiltered.fq.gz  

### Assembled Haplotypes:  
/share/pool/CompGenomVert/RefGenomes/fPhoPho  
fPhoPho.hap1.20220427.fa.gz  
fPhoPho.hap2.20220427.fa.gz  

### pacbio adapters
Pacbio adapters can be found and saved with script [] and must be renamed to `pacbio_vectors_db` and placed inside the folder `DATA/HiFiAdapterFilt/DB/pacbio_vectors_db`