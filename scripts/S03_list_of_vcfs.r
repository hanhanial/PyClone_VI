'
script to prepare input for S05 nextflow

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

idir = paste0(wdir,"/S02b_mpileup_CodingRegs/")

vcfs = list.files(idir)
vcfs = tibble(DNA_lib = vcfs,
              VCF_dir = paste0(idir,vcfs)) %>% 
  mutate(DNA_lib = gsub(pattern = ".vcf",replacement = "",x = DNA_lib)) %>% 
  mutate(coding_VCF = paste0(VCF_dir,"/mpileup_Coding.vcf"),
         noncoding_VCF = paste0(VCF_dir,"/mpileup_Noncoding.vcf")) %>% 
  inner_join(mf %>% select(Unified_ID,DNA_lib)) %>% 
  select(Unified_ID,DNA_lib,coding_VCF,noncoding_VCF)

write_csv(vcfs,paste0(wdir,"S03_list_of_vcfs.csv"))
