# **Genome Assembly Pipeline**



Table of Contents
========
- [Input](#input)
- [Dependencies](#dependencies)
- [Testing](#testing)
- [Environment Setup](#environments)
- [Running the Pipeline](#run-pipeline)
- [Pipeline Visualization](#pipeline-visualization)



## Input
- Subreads or fasta hifi reads that must be consolidated into a single file.
- HiC reads should be in two files, {filename}_R1.{fastq|fastq.gz|fq.gz} and {filename}_R2.{fastq|fastq.gz|fq.gz}.

## Dependencies

The following packages need to be set up separately. They can be installed as per the instructions on their respective GitHub pages or in a manner suitable for your system. If your setup differs significantly, you may need to manually adjust the calls to these programs within *`makefiles/STATS.smk`*:

In this pipeline we use them with modules which can be changed in the [*config.yaml*](config/config.yaml) file.

- [GATK 4.4.0.0](https://github.com/broadinstitute/gatk)
- [SALSA2 2.3](https://github.com/marbl/SALSA)


### Genome Profile

This pipeline creates a kmer based statistics using the input PacBio HiFi reads. Idealy you want PCR-free libraries for these reads to reduce bias. 
You also want high enough coverage to separate low-coverage kmers from real error kmers. 
From Genescope you would want the error peak to not overlap the heterozygous peak.

Because of the different error rates and sequencing biases of the different libraries it is generally not recommended to combine different datasets (HiFi, HiC, Illumina short reads) datasets. However, it is recommended to do separate genome profiles for each dataset.

After the pipeline runs and you have your genome profile, make sure to see how well the model fits your data and consider changing the xxx parameters to help the model fit your data better.

### HiFiasm

You can incorporate HiC reads into the Hifiasm run or not regardless of whether you have it by enabling the option INCLUDE_HIC. If you run HiFiasm with and without HIC reads make sure to change the prefix of your run ("file_prefix") so as not to overwrite the files or cause errors because of preexisting files.

### Flye

Given that Flye is easier to run and has fewer options to choose from it is implemented into a separate script found [here](flye.sh). Commented out are two self-explanatory options : --scaffold --keep-haplotypes, which can be used if wanted/needed.

There is also a snakefile for running FLYE but is unused and is given only for inspiration on how it could work.

### Statistics

In this directory [(here)](./stats) you can also find a standalone script for doing post-assembly analysis on either flye or subsequent assemblies (curated, polished, etc)

### Polishing

In here you can also find a folder with polishing scripts that use deepvariant, hypo or nextplish independently. You can then pass any of them as input to the other and use the scripts under [stats](./stats/) for analysis.

### BUSCO

For [BUSCO](https://busco.ezlab.org/), it's recommended to use it `offline`. Download your relevant lineage from [here](https://busco-data.ezlab.org/v5/data/lineages/), unpack and place in a folder:

```
curl -O <link>
tar -xf <file> # becomes <folder>
mkdir busco_downloads
mkdir busco_downloads/lineages
mv <folder> busco_downloads/lineages
```

Specify the path in the [*config.yaml*](config/config.yaml) file.

### Deepconsensus
[Deepconsensus](https://github.com/google/deepconsensus) requires a directory with a model specific to the version used (e.g., 1.2). Download the models from here and place them inside the model folder:

```
wget --no-check-certificate https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/models/v1.2/model_checkpoint/checkpoint.index  ./model/
```

You later need to provide the checkpoint folder  of this model inside the [*config.yaml*](config/config.yaml) file.

### Filter PacBio Adapters
Find and save PacBio adapters using the custom script *`utils/pbadapterfilt.sh`*. Rename the result to *`pacbio_vectors_db`* and place it inside the *`DATA/HiFiAdapterFilt/DB/`* folder.

### Salsa2

We use 2 scripts provided by Salsa2. Download Salsa locally and provide the pipeline with a path to the bin folder.

### Kraken decontamination

This pipeline also always runs a decontamination step with kraken2.
You can download a database for kraken2 from [here](https://benlangmead.github.io/aws-indexes/k2) and specify it in the [*config.yaml*](config/config.yaml) file together with a confidence score for choosing which reads will be considered classified and unclassified.

Choosing a confidence score is a non-trivial decision. Feel free to change the confidence score and use snakemake's `"--until"` option to rerun the decontamination.

## Testing

## Entire pipeline:

To test the entire pipeline you can download data from published assemblies with small genomes.
Note that the input of ccs are circular, raw pacbio reads in `.bam` format.
Example assembly raw data:

`https://www.ebi.ac.uk/ena/browser/view/PRJNA765108`

## Parts of the pipeline
### Process raw subreads
You can also test the first part of the pipeline, the creation of hifi `.fastq` reads from raw PacBio reads with the 
deepconsensus test data provided from here:

```
https://console.cloud.google.com/storage/browser/brain-genomics-public/research/deepconsensus/test_data?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false
```

For more info on how to test this part of the pipeline see the deepconsensus github:

`https://github.com/google/deepconsensus/blob/r1.2/docs/quick_start.md`

### HiC editing
 
This pipeline also includes an optional file for processing HiC reads, *`PROCESS_HIC_READS.smk`*. Test this by subsampling your own HiC paired-end Illumina reads with `seqtk sample` or by using data from a public database.

The purpose of this pipeline is to trim, decontaminate and error correct the HiC files and estimate how much they been improved. The user can then use the version of the HiC reads they prefer.

### Choose right hifiasm options for best assembly


It is recommended to begin by testing different hifiasm options before performing purge_dups and hic scaffolding (if available) and assembly visialization. More importance should be placed here in better separating the haplomes and improcing BUSCO scores, since contiguity can be improved by Scaffolding and more deduplication.


To do this, you can set a list of hifiasm options in the [*config.yaml*](config/config.yaml) file under 
HIFIASM_OPTIONS and HIFIASM_PURGE. You can then run this part of the pipeline and produce statistics for the different hifiasm options by setting the following options to True: 

DO_CCS, DO_ADAPT_FILT (optional)

DO_BUSCO, DO_STATS

We can then select the hifiasm options we want (comment them out of [*config.yaml*](config/config.yaml) file) and run the second part of the pipeline with only 1 hifiasm pair of haplomes.

## Environment(s)

In order to run the pipeline, you need a basic conda environment with python>=3.9 and mamba on your system. 
You can create that with:

```
conda create -y -n genome_assembly "python>=3.9" "mamba>=0.22.1" -c conda-forge -c bioconda -c defaults
conda activate genome_assembly
```

You can then install the environments used by the different rules of the pipeline in a terminal with an internet connection:

```
# conda activate genome_assembly
snakemake --snakefile snakefile snakefile.smk--cores 5 -p -r -w 5 --verbose --use-conda  --conda-create-envs-only
# second pipeline for processing HiC data
snakemake --snakefile PROCESS_HIC_READS.smk --cores 5 -p -r -w 5 --verbose --use-conda  --conda-create-envs-only
```

This can take several minutes and needs to be done every time the environment files are changed.

## Test run:

```
snakemake -s snakefile.smk --dry-run --cores 5 -p -r -w 5 --verbose
```

## Run pipeline:
You can then run the pipeline with:

```
snakemake \
    --snakefile snakefile.smk \
    --keep-going \
    --latency-wait 300 \
    -j ${THREADS} \
    --default-resources "tmpdir='/path/to/tmp'" \
    --verbose \
    --use-conda \
    --use-envmodules \
    --printshellcmds \
    --reason \
    --nolock \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --stats "./stats.json" 
```

## Visualize the steps of the pipeline:

You can create a represenation of the pipeline with:
```
snakemake -s snakefile.smk --dag --forceall | dot -Tpng > graph_of_jobs.png
snakemake -s snakefile.smk --dag --until HIFIASM | dot -Tpng > graph_of_jobs_hifiasm.png
snakemake -s snakefile.smk --filegraph --forceall | dot -Tpng > filegraph_all.png
snakemake -s snakefile.smk --rulegraph --forceall | dot -Tpng > rulegraph_all.png

snakemake -s PROCESS_HIC_READS.smk --rulegraph --forceall | dot -Tpng > rulegraph_hic_editing.png
```

---

### First part of pipeline: process raw subreads

![From rawreads to assembly](./graphs/rulegraph_all_ccs.png) 

In this part of the pipeline you can avoid editing the subreads, if you already have fastq HiFi reads, by setting DO_CCS to False in the [*config.yaml*](config/config.yaml) file. 

---

### HiC editing pipeline

![HiC editing and QC](./graphs/rulegraph_all_scaf.png)

By setting the following optios to True and choosing the correct amount of hifiasm options you can run the deduplication and HiC scaffolding above. All steps will always run for both haplotypes.

## General methodology for assembly:

- Run Fastk and Genescope to find heterozygous peak and assess the genome profile
- Use homozygous peak as well as other hifiasm options to find best combination for assembly
- Test this option for hifiasm that includes or not the HiC
- Use best options for purge_dups and yahs scaffolding
- Polish
- Decontaminate ASM with blobtools/FCS
- Remove possible mitochondrial contigs with MitoHiFi
- Manual curation of remaining contigs


# References of tools used
- BUSCO
    - Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654
    - Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323
- CCS
    - https://ccs.how/
- ACTC
    - https://github.com/PacificBiosciences/actc
- googe/deepconsensus
    - Baid, G., Cook, D.E., Shafin, K. et al. DeepConsensus improves the accuracy of sequences with a gap-aware sequence transformer. Nat Biotechnol 41, 232–238 (2023). https://doi.org/10.1038/s41587-022-01435-7
- HiFiAdapterFilt
    - https://github.com/sheinasim/HiFiAdapterFilt
- Kraken2
    - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
- BWA
    - Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754-1760. [PMID: 19451168]. (if you use the BWA-backtrack algorithm)
    - Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics, 26, 589-595. [PMID: 20080505]. (if you use the BWA-SW algorithm)
    - Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN]. (if you use the BWA-MEM algorithm or the fastmap command, or want to cite the whole BWA package)
- Samtools
    - Twelve years of SAMtools and BCFtools
    Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
- Salsa2
    - Ghurye, J., Pop, M., Koren, S., Bickhart, D., & Chin, C. S. (2017). Scaffolding of long read assemblies using long range contact information. BMC genomics, 18(1), 527. [Link](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3879-z)
    - Ghurye, J., Rhie, A., Walenz, B.P., Schmitt, A., Selvaraj, S., Pop, M., Phillippy, A.M. and Koren, S., 2018. Integrating Hi-C links with assembly graphs for chromosome-scale assembly. bioRxiv, p.261149 [Link](https://www.biorxiv.org/content/early/2018/02/07/261149)
- GATKv4.4
    - Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
- YaHs
    - Chenxi Zhou, Shane A McCarthy, Richard Durbin, YaHS: yet another Hi-C scaffolding tool, Bioinformatics, Volume 39, Issue 1, January 2023, btac808, https://doi.org/10.1093/bioinformatics/btac808
- Hicstuff
    - Cyril Matthey-Doret, Lyam Baudry, Amaury Bignaud, Axel Cournac, Remi-Montagne, Nadège Guiglielmoni, Théo Foutel Rodier and Vittore F. Scolari. 2020. hicstuff: Simple library/pipeline to generate and handle Hi-C data . Zenodo. http://doi.org/10.5281/zenodo.4066363
- Hifiasm
    - Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021) Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nat Methods, 18:170-175. https://doi.org/10.1038/s41592-020-01056-5
    - Cheng, H., Jarvis, E.D., Fedrigo, O., Koepfli, K.P., Urban, L., Gemmell, N.J., Li, H. (2022) Haplotype-resolved assembly of diploid genomes without parental data. Nature Biotechnology, 40:1332–1335. https://doi.org/10.1038/s41587-022-01261-x
- Purge-dups
    - Dengfeng Guan, Shane A McCarthy, Jonathan Wood, Kerstin Howe, Yadong Wang, Richard Durbin, Identifying and removing haplotypic duplication in primary genome assemblies, Bioinformatics, Volume 36, Issue 9, May 2020, Pages 2896–2898, https://doi.org/10.1093/bioinformatics/btaa025
- minimap2
    - Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191
    - Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705
- seqkit
    - Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
- FastK
    - https://github.com/thegenemyers/FASTK
- GenescopeFK, MerquryFK, KatGC
    - https://github.com/thegenemyers/GENESCOPE.FK
- GAAS
    - https://github.com/NBISweden/GAAS
- QUAST
    - Alla Mikheenko, Andrey Prjibelski, Vladislav Saveliev, Dmitry Antipov, Alexey Gurevich,
    Versatile genome assembly evaluation with QUAST-LG,
    Bioinformatics (2018) 34 (13): i142-i150. doi: 10.1093/bioinformatics/bty266
    First published online: June 27, 2018
- GenomeTools
    - G. Gremme, S. Steinbiss and S. Kurtz.
    GenomeTools: a comprehensive software library for efficient processing of structured genome annotations.
    IEEE/ACM Transactions on Computational Biology and Bioinformatics 2013, 10(3):645–656
- Gfastar
    - Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs 
    Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis doi: https://doi.org/10.1093/bioinformatics/btac460
