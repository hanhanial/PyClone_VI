
# MAF file
maf=$1

# output file
ofile=$2


# HG19 coding definitions
#echo -e "Hugo_Symbol\tchrom\tchromStart\tchromEnd\tVariant_Classification\tVariant_Type\tt_alt_count\tt_ref_count" > $ofile && cat "$maf" | grep -v "^#" | awk -F '\t' '$9~/\yDe_novo_Start_InFrame\y|\yDe_novo_Start_OutOfFrame\y|\yStart_Codon_Ins\y|\yStart_Codon_Del\y|\yStart_Codon_SNP\y|\yFrame_Shift_Del\y|\yFrame_Shift_Ins\y|\yIn_Frame_Del\y|\yIn_Frame_Ins\y|\yMissense_Mutation\y|\yNonsense_Mutation\y|\yNonstop_Mutation\y|\ySilent\y|\ySplice_Site\y|\yTranslation_Start_Site\y/{ print $1"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10"\t"$81"\t"$82 }' >> $ofile

# HG38 coding definitions
echo -e "Hugo_Symbol\tchrom\tchromStart\tchromEnd\tVariant_Classification\tVariant_Type\tt_alt_count\tt_ref_count" > $ofile && cat "$maf" | grep -v "^#" | awk -F '\t' '$9~/\DE_NOVO_START_IN_FRAME\y|\yDE_NOVO_START_OUT_FRAME\y|\ySTART_CODON_INS\y|\ySTART_CODON_DEL\y|\ySTART_CODON_SNP\y|\yFrame_Shift_Del\y|\yFrame_Shift_Ins\y|\yIn_Frame_Del\y|\yIn_Frame_Ins\y|\yMissense_Mutation\y|\yNonsense_Mutation\y|\yNonstop_Mutation\y|\ySilent\y|\ySplice_Site\y|\yTranslation_Start_Site\y/{ print $1"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10"\t"$81"\t"$82 }' >> $ofile


# 2 hot spots for TERT promoters
echo -e "TERT\tchr5\t1295113\t1295113\tNA\tSNP\tNA\tNA\nTERT\tchr5\t1295135\t1295135\tNA\tSNP\tNA\tNA" >> $ofile
