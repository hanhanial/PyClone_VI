
# MAF file
maf=$1

# output file
ofile=$2


echo -e "Hugo_Symbol\tchrom\tchromStart\tchromEnd\tVariant_Classification\tVariant_Type\tt_alt_count\tt_ref_count" > $ofile && cat "$maf" | grep -v "^#" | awk -F '\t'  '$10~/\ySNP\y/{print $1"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10"\t"$81"\t"$82}' >> $ofile

