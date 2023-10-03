# find where the adapters used in pacbio sequencing
# are mentioned and extract the dna sequence
# to check for leftover adapters in the reads 
# tags: deepconsensus, adapter contamination

DNAseqs=$(grep "Adaptor" $1 | grep -oE  ">{1}[ACGT]+")
ID=$(grep -oE "m[0-9a-z]+_[0-9a-z]+_[0-9a-z]+" <<< $(echo $1))

count=1
for seq in $DNAseqs; do 
echo "PacBio_Adapter_${count}" >> ./pacbio_adapters_${ID}
echo "${seq}" >> ./pacbio_adapters_${ID}
count=$(( ${count} + 1 ))
done
