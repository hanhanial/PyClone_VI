'
Script: S01a
get bed file of coding regions using MAF files (funcotator results)
--> from MAF files: 
for each sample, keep mutations with coding definition
(i.e. Variant_Classification in coding_def) and
mutations in TERT promoter region 
(i.e. in chromosome 5, between position 1294990 and position 1299990)
--> then use chromosome, start, end of those mutations as bed file for the sample
'

suppressMessages(library(tidyverse))
suppressMessages(library(doMC))

args = commandArgs(trailingOnly = T)

# master file
mf = args[1]

# working directory
wdir = args[2]

'
mf = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/pyclone_masterfile.csv"
wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210318/"
' 

mf = read_csv(mf)

# output directory to store DNA libs' bed files of coding regions from oncotator
odir = paste0(wdir,"/S01a_mutect_CodingRegs_perSam/")
system(paste0("mkdir -p ",odir))

mf = mf %>% 
  filter(!is.na(funcotator_dir))

script1 = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts/mutect_SNPs_perSam.sh"
script2 = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts/mutect_CodingRegs_perSam.sh"


registerDoMC(10)
final_out = foreach (d=seq(from=1,to=nrow(mf),by=1)) %dopar% {
  # d = 1
  
  dna_lib = mf$DNA_lib[d]
  
  print(paste0("Processing DNA lib ",d,": ",dna_lib))
  
  maf = mf$funcotator_dir[d]
  
  # get all SNPs
  cmd = paste0(script1," ",
               maf," ",
               odir,dna_lib,"_SNPs.bed")
  system(cmd)
  
  # get all coding mutations
  cmd = paste0(script2," ",
               maf," ",
               odir,dna_lib,"_coding.bed")
  system(cmd)
  
  
  
  return("done!")
  
}
