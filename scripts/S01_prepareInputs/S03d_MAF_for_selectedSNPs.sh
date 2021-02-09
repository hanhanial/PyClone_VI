# to get MAF for selected SNPs (i.e. SNPs input for PyClone)
selectedSNPs=$1
maf=$2
ofile=$3
grep -F -f <(awk '$1!~/^#/{print $1,$2}' OFS="\t" $selectedSNPs) $maf >$ofile
