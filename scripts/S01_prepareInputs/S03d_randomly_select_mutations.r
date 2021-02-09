'
Script: S03d
for each patients, randomly select mutations from mpileup VCF
'

suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly = T)

# master file
# mf = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/pyclone_masterfile.csv"
mf = args[1]
mf = read_csv(mf)

# working directory
# wdir = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/output/output_20201008/S01_prepareInputs/"
wdir = args[2]

# number of mutations
# num = 1500
num = args[3]

# copy neutral
# copy_neutral = "no"
copy_neutral = args[4]

# random seed for selection mutations
# random_seed = 123
random_seed = args[5]

# input directory
if (copy_neutral == "no") {
  idir = paste0(wdir,"/S03b_mpileup_CodingRegs/")
} else {
  idir = paste0(wdir,"/S03c_mpileup_CodingRegs_CopyNeutral/")
}


# output directory 
odir = paste0(wdir,"/S03d_randomly_select_mutations/")
system(paste0("mkdir -p ",odir))

# list of normal DNA lib (patient) to run
dna_libs = list.files(idir)
dna_libs = intersect(dna_libs,mf$DNA_lib)


# input directory
if (copy_neutral == "no") {
  vcfs = tibble(DNA_lib = dna_libs,
                vcf_coding = paste0(idir,dna_libs,"/mpileup_Coding.vcf"),
                vcf_noncoding = paste0(idir,dna_libs,"/mpileup_Noncoding.vcf"))
} else {
  vcfs = tibble(DNA_lib = dna_libs,
                vcf_coding = paste0(idir,dna_libs,"/Coding_CopyNeutral.vcf"),
                vcf_noncoding = paste0(idir,dna_libs,"/Noncoding_CopyNeutral.vcf"))
}


sh_code = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PLANET_Pyclone_pipeline/scripts/S01_prepareInputs/S03d_randomly_select_mutations.sh"

if (copy_neutral == "yes") {
  odir1 = paste0(odir,"copyneutral_random_",num,"mutations/")
  print(paste0("Create VCF with ",num," randomly selected COPY NEUTRAL mutations..."))
} else {
  odir1 = paste0(odir,"random_",num,"mutations/")
  print(paste0("Create VCF with ",num," randomly selected mutations..."))
}
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


