'
Script: S01a
since phyloWGS assume all sequenza results are in the same folder
--> create soft links 
and, get a list of copy-neutral segments for each DNA lib
'

suppressMessages(library(tidyverse))
suppressMessages(library(doMC))

args = commandArgs(trailingOnly = T)

# master file
# mf = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/pyclone_masterfile.csv"
mf = args[1]
mf = read_csv(mf)

# working directory
# wdir = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/output/output_20201008/S01_prepareInputs/"
wdir = args[2]

# directory to store copy-neutral regions of each sample
outp = paste0(wdir,"/S01a_sequenza_CopyNeuRegs_perSam/")
system(paste0("mkdir -p ",outp))

sequenzas = mf %>% 
  select(Unified_ID,DNA_lib,sequenza_dir,CN) %>% 
  filter(!is.na(sequenza_dir))


# for (i in 1:nrow(sequenzas))
# {
registerDoMC(32)
final_out = foreach (i=seq(from=1,to=nrow(sequenzas),by=1)) %dopar% {
  # i = 1
  
  pat = sequenzas$Unified_ID[i]
  dna_lib = sequenzas$DNA_lib[i]
  
  print(paste0("Processing sample ",i,": ",
               dna_lib))
  
  odir = paste0(outp,dna_lib)
  mkdir_cmd = paste0("mkdir -p ",odir)
  system(mkdir_cmd)
  
  # create softlinks for some files
  files_to_softlink = paste0(dna_lib,c("_confints_CP.txt","_alternative_solutions.txt"))
  files_to_softlink = paste0(sequenzas$sequenza_dir[i],"/",
                             files_to_softlink)
  cmd = paste0("ln -s -f ",paste0(files_to_softlink,collapse = " ")," -t ",odir)
  system(cmd)
  
  
  # for segment files, keep only those with copy neutral regions
  segs = read.table(paste0(sequenzas$sequenza_dir[i],
                           "/",dna_lib,"_segments.txt"),
                    header = T)
  segs = segs %>% 
    mutate(seg_length = end.pos - start.pos) %>% 
    select(chromosome,start.pos,end.pos,seg_length,everything())
  
  # CN is the minor and major allele copy number in copy-neutral regions
  CN = sequenzas$CN[i]
  
  # keep only segments in copy-neutral regions
  segs1 = segs %>% 
    filter(A==CN & B==CN)
  
  write.table(segs1,paste0(odir,"/",dna_lib,"_segments.txt"),
              quote = F,row.names = F,sep = "\t")
  return(segs1)
}
