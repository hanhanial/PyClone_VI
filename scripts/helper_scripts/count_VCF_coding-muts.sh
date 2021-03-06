dna_lib=$1
dir=$2
ofile=$3

num_coding=`cat "$dir"/"$dna_lib"/mpileup_Coding.vcf | grep -v "^#" | wc -l`
num_noncoding=`cat "$dir"/"$dna_lib"/mpileup_Noncoding.vcf | grep -v "^#" | wc -l`


if [ ! -f "$dir"/"$ofile" ]; then
  # if $dir/$ofile does not exist, create it and write its header line
  echo -e "DNA_lib\tnum_coding\tnum_noncoding" > "$dir"/"$ofile"
  echo -e ""$dna_lib"\t"$num_coding"\t"$num_noncoding"" >> "$dir"/"$ofile"
else
  echo -e ""$dna_lib"\t"$num_coding"\t"$num_noncoding"" >> "$dir"/"$ofile"
fi

