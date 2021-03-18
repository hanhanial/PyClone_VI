#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2018-07-31 17:07:02
# File Name: var_total_reads_from_vcf.py
# Description: 

# 2020-10-05 edited by Hannah, to take into account muts with no AD in mpileup
#########################################################################

import os
import sys
import argparse
import gzip

def get_mutation_form(mutectDir=None,samples=[]):
    mutations={}
    for smp in samples:
        with open('{}/{}/{}-mutect2.PASSED.vcf'.format(mutectDir,smp,smp),'rt') as input:
            for line in input:
                if line.startswith('#'):
                    pass
                else:
                    line=line.rstrip()
                    fields=line.split()
                    mutations['_'.join(fields[:2])]={'ref':fields[3].upper(),'alt':fields[4].upper()}
    return mutations

def main():
    parser=argparse.ArgumentParser(
        description='Basing the ref/alt information from mutect vcf, this '+\
        'script can count number of alt/all reads from mpileup vcf')
    parser.add_argument('-v','--vcf',required=True,
        help='mpileup vcf file for multiple tumor sectors')
    parser.add_argument('-p','--patient',required=True,
        help='the patient id to be used as the prefix of output files')
    parser.add_argument('-m','--mutectDir',required=True,
        help='the base directory of mutect results')

    args=parser.parse_args()
    
    with open(args.vcf,'r') as vcf, \
        open(args.patient+'.sectors','w') as sectors_f, \
        open(args.patient+'_mutant_reads.tsv','w') as alt_tsv, \
        open(args.patient+'_read_depth.tsv','w') as cov_tsv:
        sectors=[]
        mutations={}
        for line in vcf:
            line=line.rstrip()
            if line.startswith('##'):
                pass
            elif line.startswith('#CHROM'):
                sectors=line.split('\t')[10:]
                sectors=[x[-6:] for x in sectors]
                sectors_f.write('\n'.join(sectors)+'\n')
                alt_tsv.write('Chromosome\tPosition\tChange\tGene\t'+'\t'.join(sectors)+'\n')
                cov_tsv.write('Chromosome\tPosition\tChange\tGene\t'+'\t'.join(sectors)+'\n')
                mutations=get_mutation_form(mutectDir=args.mutectDir,samples=sectors)
            else:
                fields=line.split('\t')
                chr=fields[0]
                pos=fields[1]
                ref=fields[3]
                alt=fields[4]
                if alt=='.':
                    continue
                alleles=[ref]+[x for x in alt.split(',') if len(x)==1]
                alleles=[x.upper() for x in alleles]

                try:
                    mutect_ref=mutations['_'.join([chr,pos])]['ref']
                    mutect_alt=mutations['_'.join([chr,pos])]['alt']
                except KeyError:
                    continue
                fmt=fields[8]
                ref_ADs=[]
                alt_ADs=[]
                total_cov=[]
                AD_index=fmt.split(':').index('AD')
                try:
                    ref_index=alleles.index(mutect_ref)
                    alt_index=alleles.index(mutect_alt)
                except ValueError:
                    print('In mutect, REF is {} and ALT is {}. But at least one is missing in the mpileup VCF:'.format(mutect_ref,mutect_alt)) 
                    print(line+'\n')
                    continue
                
                for i in range(0,len(sectors)):
                    smp_info=fields[10+i]
                    try:
                        ref_ADs.append(smp_info.split(':')[AD_index].split(',')[ref_index])
                        alt_ADs.append(smp_info.split(':')[AD_index].split(',')[alt_index])
                        total_cov.append(int(ref_ADs[-1])+int(alt_ADs[-1]))
                    except IndexError:
                        print('This mutation doesnt have allelic depths: {}'.format(smp_info))        
                        continue
                alt_tsv.write('\t'.join([chr,pos,ref+'>'+alt,'Unknown']+alt_ADs)+'\n')
                cov_tsv.write('\t'.join([chr,pos,ref+'>'+alt,'Unknown']+[str(x) for x in total_cov])+'\n')

if __name__=='__main__':
    main()


