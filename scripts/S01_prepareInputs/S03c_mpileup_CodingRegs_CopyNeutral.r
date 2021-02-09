'
Script: S03c
for each DNA lib, 
using copy-neutral bed files (from sequenza results) and 
mpileup of coding mutations and 
mpileup of non-coding mutations

--> output 2 files: 
mutations in copy-neutral + coding regions
mutations in copy-neutral + non-coding regions
'

suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly = T)

# master file
# mf = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/pyclone_masterfile.csv"
mf = args[1]
mf = read_csv(mf)
mf = mf %>% 
  filter(Sample_type=="Normal") %>% 
  select(Unified_ID,DNA_lib)

# working directory
# wdir = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/output/output_20201008/S01_prepareInputs/"
wdir = args[2]

# directory to mpileup's mutations in coding regions (i.e. output directory of S03b.)
vcf_dir = paste0(wdir,"/S03b_mpileup_CodingRegs/")
dna_libs = list.files(vcf_dir)
dna_libs = tibble(DNA_lib = dna_libs,
                  vcf_dir = paste0(vcf_dir,dna_libs))


# directory to DNA libs' bed files of copy-neutral regions (i.e. output directory of S01b.)
bed_dir = paste0(wdir,"/S01b_sequenza_CopyNeuRegs_perPat/")
bed_files = list.files(bed_dir)
bed_files = tibble(Unified_ID = bed_files,
                   bed_dir = paste0(bed_dir,bed_files)) %>% 
  mutate(Unified_ID = gsub(pattern = "_CopyNeutralSegments.bed",replacement = "",x = Unified_ID)) %>% 
  inner_join(mf)

dna_libs = dna_libs %>% 
  inner_join(bed_files)

outp = paste0(wdir,"/S03c_mpileup_CodingRegs_CopyNeutral/")
system(paste0("mkdir -p ",outp))

num_muts = tibble(DNA_lib = character(),
                  num_coding_copyneutral = numeric(),
                  num_noncoding_copyneutral = numeric())
# remove existing number_of_CopyNeutral-mutations.tsv file that counts number of muts
system(paste0("rm -f ",outp,"/number_of_CopyNeutral-mutations.tsv"))


for (d in 1:nrow(dna_libs))
{
  # d = 1
  dna_lib = dna_libs$DNA_lib[d]
  
  print(paste0("Processing DNA lib ",d,": ",dna_lib))
  
  # create a folder for each DNA lib
  odir = paste0(outp,dna_lib)
  system(paste0("mkdir -p ",odir))
  
  copyneutral_bed = dna_libs$bed_dir[d]
  vcf_coding = paste0(dna_libs$vcf_dir[d],"/mpileup_Coding.vcf.gz")
  vcf_noncoding = paste0(dna_libs$vcf_dir[d],"/mpileup_Noncoding.vcf.gz")
  
  ##### mutations in copy-neutral + coding regions
  cmd = paste0("vcftools --gzvcf ",vcf_coding,
               " --recode --recode-INFO-all --bed ",copyneutral_bed,
               " --out ",odir,"/coding_copyneutral")
  system(cmd)
  
  # rename file
  system(paste0("mv ",odir,"/coding_copyneutral.recode.vcf ",odir,"/Coding_CopyNeutral.vcf"))
  
  # compress
  system(paste0("bgzip -c ",odir,"/Coding_CopyNeutral.vcf > ",odir,"/Coding_CopyNeutral.vcf.gz"))
  
  # index compressed file
  system(paste0("bcftools index ",odir,"/Coding_CopyNeutral.vcf.gz"))
  
  ##### mutations in copy-neutral + non-coding regions
  cmd = paste0("vcftools --gzvcf ",vcf_noncoding,
               " --recode --recode-INFO-all --bed ",copyneutral_bed,
               " --out ",odir,"/noncoding_copyneutral")
  system(cmd)
  
  # rename file
  system(paste0("mv ",odir,"/noncoding_copyneutral.recode.vcf ",odir,"/Noncoding_CopyNeutral.vcf"))
  
  # compress
  system(paste0("bgzip -c ",odir,"/Noncoding_CopyNeutral.vcf > ",odir,"/Noncoding_CopyNeutral.vcf.gz"))
  
  # index compressed file
  system(paste0("bcftools index ",odir,"/Noncoding_CopyNeutral.vcf.gz"))
  
  # number of coding muts
  cmd = paste0("/mnt/projects/lailhh/workspace/pipelines/PyClone/PLANET_Pyclone_pipeline/scripts/S01_prepareInputs/count_VCF_copyneutral-muts.sh ",
               dna_lib," ",outp," number_of_CopyNeutral-mutations.tsv")
  system(cmd)
}
