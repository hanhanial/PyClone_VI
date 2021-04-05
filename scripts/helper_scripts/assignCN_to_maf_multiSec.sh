pat_id=$1
segs_file=$2
coding_var_reads=$3
coding_depth=$4
noncoding_var_reads=$5
noncoding_depth=$6
script_dir=$7

python ${script_dir}/assignCN_to_maf_multiSec.py -m $coding_var_reads -s $segs_file -p $pat_id -d $coding_depth
cd "$pat_id"_pyclone_input
for i in *; do mv "$i" coding_"$i"; done
cd ../
  
python ${script_dir}/assignCN_to_maf_multiSec.py -m $noncoding_var_reads -s $segs_file -p $pat_id -d $noncoding_depth
cd "$pat_id"_pyclone_input
noncoding_files=$(ls | grep -v coding)
IFS=$'\n';for i in ${noncoding_files}; do mv "$i" noncoding_"$i"; done
cd ../
