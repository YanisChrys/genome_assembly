#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q large.q,medium.q
#$ -N process_hic
#$ -pe smp 21
#$ -e /share/pool/CompGenomVert/EpiAus/workflow/logfiles
#$ -o /share/pool/CompGenomVert/EpiAus/workflow/logfiles
#$ -M i.chrysostomakis@leibniz-lib.de
#$ -m beas

module load anaconda3/2022.05
conda activate genofish

export THREADS=$(expr ${NSLOTS} - 1)

snakemake \
    --snakefile makefiles/PROCESS_HIC_READS.smk \
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
    --rerun-triggers mtime \
    --rerun-incomplete \
    --stats "./stats/stats_hic_process.json" 

    #--report "./report.html"

# snakemake -s snakefile.smk --dag --forceall | dot -Tpng > graph_of_jobs.png && \
# snakemake -s snakefile.smk --filegraph --forceall | dot -Tpdf > filegraph_all.pdf  && \
# snakemake -s snakefile.smk --rulegraph --forceall | dot -Tpdf > rulegraph_all.pdf 


# snakemake -s snakefile.smk --dag --forceall | dot -Tpng > graph_of_jobs_assembly.png && \
# snakemake -s snakefile.smk --filegraph --forceall | dot -Tpdf > filegraph_all_ccs.pdf  && \
# snakemake -s snakefile.smk --rulegraph --forceall | dot -Tpng > rulegraph_all_pol.png 


