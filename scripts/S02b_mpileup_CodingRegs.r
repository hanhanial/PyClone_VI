'
Script: S02b
from mpileup VCF and coding bed file (from mutect MAF)
use vcftools to separate mutations into 2 groups:
- mutations in coding regions 
- mutations in non-oding regions 

'

suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly = T)

# master file
# mf = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/pyclone_masterfile.csv"
mf = args[1]
mf = read_csv(mf)

# working directory
# wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210318/"
wdir = args[2]

# input directory
idir = paste0(wdir,"/S02a_compress_index_mpileup_VCF/")

# directory to bed files (i.e. lists of merged coding regions) of all patients
bed_dir = paste0(wdir,"/S01b_mutect_CodingRegs_perPat/")

# output directory to store mpileup's mutations in coding/non-coding regions
outp = paste0(wdir,"/S02b_mpileup_CodingRegs/")
system(paste0("mkdir -p ",outp))

mpileups = list.files(idir)
mpileups = tibble(vcf = paste0(idir,"/",mpileups,"/SNPs.vcf.gz"),
                  DNA_lib = mpileups)
mpileups = mpileups %>% 
  inner_join(mf)

num_muts = tibble(DNA_lib = character(),
                  num_codingSNP = numeric(),
                  num_noncodingSNP = numeric())

# remove existing number_of_Coding-mutations.tsv file that counts number of muts
system(paste0("rm -f ",outp,"/number_of_Coding-mutations.tsv"))

for (i in 1:nrow(mpileups))
{
  # i = 1
  
  print(paste0("Processing DNA lib ",i))
  dna_lib = mpileups$DNA_lib[i]
  
  odir = paste0(outp,"/",dna_lib)
  system(paste0("mkdir -p ",odir))
  
  # bed file
  bed_file = paste0(bed_dir,mpileups$Unified_ID[i],"_CodingSegments.bed")
  
  ###### SNP in coding regions (from mutect MAF file)
  # keep mutations in coding regions
  vcftools_cmd = paste0("vcftools --gzvcf ",idir,"/",dna_lib,"/SNPs.vcf.gz",
                        " --recode --recode-INFO-all --bed ",bed_file,
                        " --out ",odir,"/coding")
  system(vcftools_cmd)
  
  # rename
  system(paste0("mv ",odir,"/coding.recode.vcf ",odir,"/mpileup_Coding.vcf"))
  
  # compress
  system(paste0("bgzip -c ",odir,"/mpileup_Coding.vcf > ",odir,"/mpileup_Coding.vcf.gz"))
  
  # index files
  system(paste0("bcftools index ",odir,"/mpileup_Coding.vcf.gz"))

    
  ###### SNPs in non-coding regions 
  # keep mutations in coding regions
  vcftools_cmd = paste0("vcftools --gzvcf ",idir,"/",dna_lib,"/SNPs.vcf.gz",
                        " --recode --recode-INFO-all --exclude-bed ",bed_file,
                        " --out ",odir,"/noncoding")
  system(vcftools_cmd)
  
  # rename
  system(paste0("mv ",odir,"/noncoding.recode.vcf ",odir,"/mpileup_Noncoding.vcf"))
  
  # compress
  system(paste0("bgzip -c ",odir,"/mpileup_Noncoding.vcf > ",odir,"/mpileup_Noncoding.vcf.gz"))
  
  # index files
  system(paste0("bcftools index ",odir,"/mpileup_Noncoding.vcf.gz"))
  
  ####### number of coding/noncoding muts
  cmd = paste0("/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts/count_VCF_coding-muts.sh ",
               dna_lib," ",outp," number_of_Coding-SNP.tsv")
  system(cmd)
  
}

