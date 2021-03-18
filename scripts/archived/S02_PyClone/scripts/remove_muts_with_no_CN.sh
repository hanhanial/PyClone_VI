# to remove mutations in segments with NA copy numbers (i.e. NA in A or B columns in sequenza segment output)
# --> remove all mutations with "" values in minor_cn and major_cn columns


# patient ID
pat_id=$1
pat_input_dir=${pat_id}_pyclone_input

echo $(ls "$pat_input_dir"/*.tsv)
echo $(pwd)

for sec in $(ls "$pat_input_dir"/*.tsv)
do
    ofile=tmp.tsv
    awk '($5 != "") && ($6 != "")' $sec > $ofile && mv $ofile $sec
done
