export PATH="/home/chuakp/anaconda3/bin:$PATH"

export QT_QPA_PLATFORM='offscreen'

prior=major_copy_number
iters=20000
burnin=2000
# pycEnv=/home/hfyteng/.conda/envs/pyclone/

source activate PyClone

# Sample name
sam=$1
sam_input_dir=${sam}_pyclone_input

# Purity file, make sure second column is purity
pur=$2

echo $(pwd)

mkdir -p output/${sam}

rm -f output/${sam}/purity.txt

#for sec in $(ls ${sam}*.tsv)
for sec in $(ls *.tsv)
do
	# Remove prefix
	filename=$(basename ${sec})
	# Remove extension
	secName=${filename%.*}
	echo ${secName}
	purity=$(grep -P "${secName}\t" ${pur})
	purity=($purity)
	echo ${purity[1]} >> output/${sam}/purity.txt
done

## Update 2018-7-10: Set default to binomial data as we have low coverage
#PyClone run_analysis_pipeline --density pyclone_binomial --burnin ${burnin} --prior ${prior} --num_iters ${iters} --tumour_contents $(cat output/${sam}/purity.txt) --in_files $(ls *.tsv) --working_dir output/${sam} --min_cluster_size 2 &> output/${sam}/setup_analysis.log 

## Update 2020-09-19: semi-manually run the steps of PyClone since it fails at plotting step
PyClone setup_analysis --density pyclone_binomial --prior ${prior} --num_iters ${iters} --tumour_contents $(cat output/${sam}/purity.txt) --in_files $(ls *.tsv) --working_dir output/${sam} &> output/${sam}/setup_analysis.log

PyClone run_analysis --config_file output/${sam}/config.yaml &> output/${sam}/run_analysis.log

mkdir -p output/${sam}/tables

PyClone build_table --burnin ${burnin} --config_file output/${sam}/config.yaml --out_file output/${sam}/tables/loci.tsv --table_type loci &> output/${sam}/build_table.log

PyClone build_table --burnin ${burnin} --config_file output/${sam}/config.yaml --out_file output/${sam}/tables/cluster.tsv --table_type cluster &>> output/${sam}/build_table.log

# PyClone build_mutations_file --in_file $(ls ${sam}*.tsv) --out_file output/${sam}/yaml --prior ${prior}
# 
# PyClone run_analysis --config_file output/${sample[$index_num]}/config.yaml &> output/${sample[${index_num}]}/run_analysis.log
# 
# PyClone build_table --burnin ${burnin} --config_file output/${sample[$index_num]}/config.yaml --out_file output/tables/${sample[$index_num]}/loci.tsv --table_type loci &> output/${sample[${index_num}]}/build_loci.log
# 
# PyClone build_table --burnin ${burnin} --config_file output/${sample[$index_num]}/config.yaml --out_file output/tables/${sample[$index_num]}/cluster.tsv --table_type cluster &>> output/${sample[${index_num}]}/build_loci.log
# 
# PyClone plot_loci --burnin ${burnin} --config_file output/${sample[$index_num]}/config.yaml --plot_file output/${sample[$index_num]}/similarity_matrix.pdf --plot_type similarity_matrix &> output/${sample[${index_num}]}/plot_loci.log
# 
# PyClone plot_clusters --burnin ${burnin} --config_file output/${sample[$index_num]}/config.yaml --plot_file output/${sample[$index_num]}/density.pdf --plot_type density &>> output/${sample[${index_num}]}/plot_loci.log
