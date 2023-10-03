




mkdir tmp
ls *.fastq.gz > FILES
# kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 10000x
kmc -k21 -t16 -m64 -ci1 -cs10000 @FILES kmcdb tmp
kmc_tools transform kmcdb histogram kmcdb_k21.hist -cx10000

# estimate kmer coverage thresholds
L=$(smudgeplot.py cutoff kmcdb_k21.hist L)
U=$(smudgeplot.py cutoff kmcdb_k21.hist U)
echo $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000

kmc_tools transform kmcdb -ci"$L" -cx"$U" reduce kmcdb_L"$L"_U"$U"
smudge_pairs kmcdb_L"$L"_U"$U" kmcdb_L"$L"_U"$U"_coverages.tsv kmcdb_L"$L"_U"$U"_pairs.tsv > kmcdb_L"$L"_U"$U"_familysizes.tsv


kmc_tools transform kmcdb -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv

