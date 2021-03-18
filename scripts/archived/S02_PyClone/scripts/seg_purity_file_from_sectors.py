#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2018-08-02 14:17:57
# File Name: seg_purity_file_from_sectors.py
# Description: 
#########################################################################

import os
import sys
import argparse

def main():
    parse=argparse.ArgumentParser(description='')
    parse.add_argument('-s','--sectors',required=True,metavar='FILE',
        help='sector file of a patient')
    parse.add_argument('-p','--patient',required=True,metavar='STR',
        help='the patient id to be used as the prefix of output file')
    parse.add_argument('-d','--seqenza_dir',required=True,metavar='STR',
        help='the basedir of sequenza results')
    args=parse.parse_args()

    with open(args.sectors) as input, open(args.patient+'.seg','w') as seg_file, open(args.patient+'.purity','w') as purity_file:
        
        sectors=[]
        for line in input:
            line=line.rstrip()
            fields=line.split()
            sectors.extend(fields)

        seg_file.write("Tumor_Sample_Barcode\tchromosome\tstart.pos\tend.pos\tBf\tN.BAF\tsd.BAF\tdepth.ratio\tN.ratio\tsd.ratio\tCNt\tA\tB\tLPP\n")
        purity_file.write("Tumor_Sample_Barcode\tpurity\tploidy_model\tploidy\n")
        for sec in sectors:
            with open('{}/{}/{}_segments.txt'.format(args.seqenza_dir,sec,sec)) as sec_seg:
                for line in sec_seg:
                    if not line.startswith('chromosome'):
                        seg_file.write(sec+"\t"+line.replace('"',''))

            mean_ploidy=None
            with open('{}/{}/{}_confints_CP.txt'.format(args.seqenza_dir,sec,sec)) as confints:
                nr=0
                for line in confints:
                    nr+=1
                    if nr==3:
                        line=line.rstrip()
                        mean_ploidy=float(line.split()[2])

            with open('{}/{}/{}_alternative_solutions.txt'.format(args.seqenza_dir,sec,sec)) as alternative:
                all_solutions=[]
                for line in alternative:
                    if not line.startswith('cellularity'):
                        line=line.rstrip()
                        fields=line.split()
                        fields.append(abs(float(fields[1])-mean_ploidy))
                        all_solutions.append(fields)

            picked=sorted(all_solutions,key=lambda solution: solution[-1])[0]
            purity_file.write('\t'.join([sec,picked[0],picked[1],str(mean_ploidy)])+'\n')


if __name__=='__main__':
    main()
