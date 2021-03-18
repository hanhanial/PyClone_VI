'
script to prepare input for S05 nextflow

'

suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly = T)

# master file
# mf = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/pyclone_masterfile.csv"
mf = args[1]
mf = read_csv(mf)

# wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210318/"
wdir = args[2]

# num_muts = 3000
num_muts = args[3]

idir = paste0(wdir,"/S03_randomly_select_mutations/random_",num_muts,"mutations/")


vcfs = list.files(idir)
vcfs = tibble(DNA_lib = vcfs,
              VCF = paste0(idir,vcfs)) %>% 
  mutate(DNA_lib = gsub(pattern = ".vcf",replacement = "",x = DNA_lib)) %>% 
  inner_join(mf %>% select(Unified_ID,DNA_lib))

write_csv(vcfs,paste0(wdir,"S04_list_of_vcfs.csv"))

# script = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts/var_total_reads_from_vcf.py"
# 
# for (n in 1:nrow(vcfs)) {
#   # n = 1
#   vcf_file = as.character(vcfs$VCF[n])
#   pat = as.character(vcfs$Unified_ID[n])
#   cmd = paste0("python ",script,
#                " -v ",vcf_file,
#                " -p ",pat,
#                " -m ",mutect_dir)
#   system(cmd)
# }
# 
