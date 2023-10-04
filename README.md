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

In this pipeline we use them with modules which can be changed in the *config.yaml* file.

- [FastK](https://github.com/thegenemyers/FASTK)
- [MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)
- [GENESCOPE.FK](https://github.com/thegenemyers/GENESCOPE.FK)
- [GATK 4.4.0.0](https://github.com/broadinstitute/gatk)
- [SALSA2 2.3](https://github.com/marbl/SALSA)

Especially for Merqury, in our system, it runs best when installed locally in your home directory. Initialize it with make after downloading from [here](https://github.com/thegenemyers/MERQURY.FK).

### BUSCO

For [BUSCO](https://busco.ezlab.org/), it's recommended to use it `offline`. Download your relevant lineage from [here](https://busco-data.ezlab.org/v5/data/lineages/), unpack and place in a folder:

```
curl -O <link>
tar -xf <file> # becomes <folder>
mkdir busco_downloads
mkdir busco_downloads/lineages
mv <folder> busco_downloads/lineages
```

Specify the path in the *config.yaml* file.

### Deepconsensus
[Deepconsensus](https://github.com/google/deepconsensus) requires a directory with a model specific to the version used (e.g., 1.2). Download the models from here and place them inside the model folder:

```
wget --no-check-certificate https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/models/v1.2/model_checkpoint/checkpoint.index  ./model/
```

You later need to provide the checkpoint folder  of this model inside the *config.yaml* file.

### Filter PacBio Adapters
Find and save PacBio adapters using the custom script *`utils/pbadapterfilt.sh`*. Rename the result to *`pacbio_vectors_db`* and place it inside the *`DATA/HiFiAdapterFilt/DB/`* folder.

### Kraken decontamination

This pipeline also always runs a decontamination step with kraken2.
You can download a database for kraken2 from [here](https://benlangmead.github.io/aws-indexes/k2) and specify it in the *config.yaml* file together with a confidence score for choosing which reads will be considered classified and unclassified.

Choosing a confidence score is a non-trivial decision. Feel free to change the confidence score and use snakemake's `"--until"` option to rerun the decontamination.

## Testing

## Entire pipeline:

To test the entire pipeline you can download data from published assemblies with small genomes.
Note that the input of ccs are circular, raw pacbio reads in `.bam` format.
Example assembly raw data:

`https://www.ebi.ac.uk/ena/browser/view/PRJNA765108`

## Parts of the pipeline
### First part
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
    --stats "./stats.json" \
    --report "./report.html"
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

### First part of pipeline:

![From rawreads to assembly](./graph_of_jobs_hifiasm.png) 

---

### HiC editing pipeline

![HiC editing and QC](./rulegraph_hic_editing.png)


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
