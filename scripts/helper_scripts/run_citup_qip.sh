pat=$1

clustprev=${pat}_cluster_prev.txt

maxnode=$(wc -l ${clustprev})
maxnode=($maxnode)

if [ ${maxnode[0]} -gt 10 ]
then
	maxnode=10
else
	maxnode=${maxnode[0]}
fi

patname=${ifile::4}
odir=${pat}_tmpdir
resfile=${pat}_results.h5

for i in $(seq 0 $(expr ${maxnode} - 1))
do
	echo ${i} >> ${pat}_clust_buf.txt
done

run_citup_qip.py ${clustprev} ${pat}_clust_buf.txt ${resfile} --submit local --min_nodes 1 --max_nodes ${maxnode} --nocleanup --tmpdir ${odir} --maxjobs 24
