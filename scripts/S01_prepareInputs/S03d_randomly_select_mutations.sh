#!/bin/bash

# to run interactively: ./S03d_randomly_select_mutations.sh
# to submit: qsub S03d_randomly_select_mutations.sh

# Modify according to need
#$ -S /bin/bash 
#$ -M lailhh@gis.a-star.edu.sg 
#$ -V
#$ -cwd

# Set number of CPU here 
#$ -pe OpenMP 4

#$ -l mem_free=10G,h_rt=12:00:00


# to randomly select $num_muts mutations from $vcf file, and output it to $ofile
vcf_coding=$1
vcf_noncoding=$2
num_muts=$3
ofile=$4
random_seed=$5

# function to set random seed
get_seeded_random()
{
 seed="$1"
 openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
	 </dev/zero 2>/dev/null
}


# remove existing $ofile, if there's any
rm -f "$ofile"

# get total number of coding and non-coding mutations
num_coding_muts=`cat $vcf_coding | grep -v "^#" | wc -l`
num_noncoding_muts=`cat $vcf_noncoding | grep -v "^#" | wc -l`

if [[ num_coding_muts -eq 0 && num_noncoding_muts -eq 0 ]]; then
  echo No mutations for this patient!
  exit
fi

# randomly select mutations
# if $num_coding_muts > $num_muts, then randomly select from $vcf_coding file
# if $num_coding_muts ==  $num_muts, then just return $vcf_coding file
# if $num_coding_muts < $num_muts, set num_noncoding_2bselected_muts=$num_muts-$num_coding_muts
## if $num_noncoding_2bselected_muts >= $num_noncoding_muts, then use all noncoding muts to top up to get $num_muts
## if $num_noncoding_2bselected_muts < $num_noncoding_muts, then randomly select $num_noncoding_2bselected_muts from noncoding muts to top up
if [ $num_coding_muts -gt $num_muts ]
then
  echo $num_muts mutations are randomly selected from coding mutations
  cat "$vcf_coding" | grep "^#" >> "$ofile" && cat "$vcf_coding" | grep -v "^#" | shuf -n "$num_muts" --random-source=<(get_seeded_random "$random_seed") >> "$ofile"
elif [ $num_coding_muts -eq $num_muts ]
then
  echo $num_muts mutations are all of coding mutations
  ln -s  "$vcf_coding" "$ofile"
else 
  num_noncoding_2bselected_muts="$(($num_muts-$num_coding_muts))"
  if [ $num_noncoding_2bselected_muts -lt $num_noncoding_muts ]
  then
    echo $num_coding_muts are all of coding mutations, $num_noncoding_2bselected_muts are randomly selected from noncoding mutations 
    cat "$vcf_coding" >> "$ofile" && cat "$vcf_noncoding" | grep -v "^#" | shuf -n "$num_noncoding_2bselected_muts" --random-source=<(get_seeded_random "$random_seed") >> "$ofile"
  else
    echo Selected mutations consist of all "("$num_coding_muts")" coding and all "("$num_noncoding_muts")" noncoding mutations
    cat "$vcf_coding" >> "$ofile" && cat "$vcf_noncoding" | grep -v "^#" >> "$ofile"
  fi
fi
