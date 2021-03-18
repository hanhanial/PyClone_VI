export PATH="/home/lailhh/miniconda3/bin:$PATH"
export QT_QPA_PLATFORM='offscreen'
source activate pyclone-vi

pat_id=$1
inp=$2
wdir=$3

mkdir -p "$pat_id"

num_clusters=40
num_restarts=100
density="beta-binomial"


# run pyclone-vi
pyclone-vi fit -i "$inp" -o "$pat_id"/"$pat_id".h5 -c "$num_clusters" -d "$density" -r "$num_restarts"

# output the final results from the best random restart
pyclone-vi write-results-file -i "$pat_id"/"$pat_id".h5 -o "$pat_id"/"$pat_id".tsv

