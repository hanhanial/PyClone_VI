'
Script: S02a
(run this in R cluster)
remove mutations with invalid chromosomes (i.e. those containing "_" in chromosome names) in mpileup VCF,
sort the VCF file,
compress the sorted file, then
index the compressed file
'

suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly = T)

# master file
mf = args[1]

# working directory
wdir = args[2]

'
mf = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/pyclone_masterfile.csv"
wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210429/"
'

mf = read_csv(mf)

mpileups = mf %>% 
  filter(!is.na(mpileup_dir) & Sample_type=="Normal")


# input directory - to keep only SNPs
idir = paste0(wdir,"/S01b_mutect_CodingRegs_perPat/")
  
# output directory
outp = paste0(wdir,"/S02a_compress_index_mpileup_VCF/")
system(paste0("mkdir -p ",outp))

for (i in 1:nrow(mpileups)) 
{
  # i = 1
  dna_lib = mpileups$DNA_lib[i]
  pat = mpileups$Unified_ID[i]
  mpileup_file = mpileups$mpileup_dir[i]
  
  print(paste0("Processing file ",i,": ",dna_lib))
  
  odir = paste0(outp,"/",dna_lib,"/")
  system(paste0("mkdir -p ",odir))
  
  #### remove mutations with invalid chromosomes (i.e. those containing "_" in chromosome names)
  cmd = paste0("/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts/remove_invalid_chr.sh ",
               mpileup_file," ",odir,"/",dna_lib,"_cleaned_muts.vcf")
  system(cmd)
  
  
  #### remove samples with low purity from mpileup VCF file
  pat_tum = mf %>% 
    filter(Unified_ID==pat & Sample_type!="Normal") %>% 
    pull(DNA_lib)
  
  cmd = paste0("bcftools view ",
               odir,"/",dna_lib,"_cleaned_muts.vcf", 
               " -s ", paste(c(dna_lib,pat_tum),collapse = ","),
               " -o ", odir,"/",dna_lib,"_cleaned.vcf")
  system(cmd)
  
  # remove the _cleaned_muts.vcf file
  system(paste0("rm ",odir,"/",dna_lib,"_cleaned_muts.vcf"))
  
  
  #### sort mpileup file
  system(paste0("cat ", odir,"/",dna_lib,"_cleaned.vcf", 
                "| vcf-sort > ",odir,dna_lib,".vcf"))
  
  # remove the unsorted file
  system(paste0("rm ",odir,"/",dna_lib,"_cleaned.vcf"))
  
  
  #### compress the sorted file
  system(paste0("bgzip -c ",odir,dna_lib,".vcf > ",
                odir,dna_lib,".vcf.gz"))
  
  
  #### index the compressed file
  system(paste0("bcftools index ",odir,dna_lib,".vcf.gz"))
  
  
  
  #### keep only SNPs
  snps_bed = paste0(idir,"/",pat,"_SNPSegments.bed")
  
  # keep SNP mutations 
  vcftools_cmd = paste0("vcftools --gzvcf ",odir,"/",dna_lib,".vcf.gz",
                        " --recode --recode-INFO-all --bed ",snps_bed,
                        " --out ",odir,"/SNPs")
  system(vcftools_cmd)
  
  # rename
  system(paste0("mv ",odir,"/SNPs.recode.vcf ",odir,"/SNPs.vcf"))
  
  # compress
  system(paste0("bgzip -c ",odir,"/SNPs.vcf > ",odir,"/SNPs.vcf.gz"))
  
  # index files
  system(paste0("bcftools index ",odir,"/SNPs.vcf.gz"))
  
  
}