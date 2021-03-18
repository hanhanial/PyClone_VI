'
Script: S03
for each patients, randomly select mutations from mpileup VCF
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

# number of mutations
# num = 3000
num = args[3]

# random seed for selection mutations
# random_seed = 123
random_seed = args[4]

# input directory
idir = paste0(wdir,"/S02b_mpileup_CodingRegs/")


# output directory 
odir = paste0(wdir,"/S03_randomly_select_mutations/")
system(paste0("mkdir -p ",odir))

# list of normal DNA lib (patient) to run
dna_libs = list.files(idir)
dna_libs = intersect(dna_libs,mf$DNA_lib)


# input directory
vcfs = tibble(DNA_lib = dna_libs,
              vcf_coding = paste0(idir,dna_libs,"/mpileup_Coding.vcf"),
              vcf_noncoding = paste0(idir,dna_libs,"/mpileup_Noncoding.vcf"))


sh_code = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts/randomly_select_mutations.sh"

odir1 = paste0(odir,"random_",num,"mutations/")
print(paste0("Create VCF with ",num," randomly selected mutations..."))
system(paste0("mkdir -p ",odir1))

for (v in 1:nrow(vcfs))
{
  # v = 1
  dna_lib = vcfs$DNA_lib[v]
  
  print(paste0("Processing file ",v,": ",dna_lib))
  
  
  cmd = paste0(sh_code," ",
               vcfs$vcf_coding[v]," ",
               vcfs$vcf_noncoding[v]," ",
               num," ",
               paste0(odir1,dna_lib,".vcf"))
  system(cmd)
  
  
}


