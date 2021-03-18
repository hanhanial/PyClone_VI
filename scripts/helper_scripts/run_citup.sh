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

run_citup_iter.py ${clustprev} ${resfile} --submit local --min_nodes 1 --max_nodes ${maxnode} --nocleanup --tmpdir ${odir} --maxjobs 24
