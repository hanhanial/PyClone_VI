'
Script: S01b
for each patient,
intersect copy-neutral segments from their samples 
--> those interesteced regions will be used as bed file to filter mutations in mpileup
'

library(GenomicRanges)
library(tidyverse)

args = commandArgs(trailingOnly = T)

# master file, or list of DNA libs to be processed for each patient
# mf = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/pyclone_masterfile.csv"
mf = args[1]
mf = read_csv(mf)

# working directory
# wdir = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S03_PyClone/output/output_20201008/S01_prepareInputs/"
wdir = args[2]

# directory storing processed sequenza outputs (i.e. copy-neutral regions) 
idir = paste0(wdir,"/S01a_sequenza_CopyNeuRegs_perSam/")

# output directory to store copy neutral regions of each patient (i.e. by intersecting copy-neutral regions of all samples of the patient) 
odir = paste0(wdir,"/S01b_sequenza_CopyNeuRegs_perPat/")
system(paste0("mkdir -p ",odir))

# DNA libs and directory to the file that storing their copy-neutral segments
copy_neutral_segs = list.files(idir)
copy_neutral_segs = tibble(DNA_lib = copy_neutral_segs,
                           seg_dir = paste0(idir,copy_neutral_segs)) %>% 
  mutate(seg_dir = paste0(seg_dir,"/",DNA_lib,"_segments.txt"))

copy_neutral_segs = copy_neutral_segs %>% 
  filter(DNA_lib %in% mf$DNA_lib) %>% 
  inner_join(mf)

pats = mf %>% 
  pull(Unified_ID) %>% 
  unique()

# to record number of (intersected) copy-neutral regions for each patient
num_copyneutral_segs = tibble(Unified_ID = character(),
                              num_copyneutral_segs = numeric())

for (p in 1:length(pats))
{
  # p = 1
  pat = pats[p]
  
  print(paste0("Processing patient ",p,": ",pat))
  
  pat_segs = copy_neutral_segs %>% 
    filter(Unified_ID==pat) %>%
    pull(seg_dir)
  
  # list of segment data of all the tumor samples of the patient
  segs = lapply(pat_segs, function(x) {dat = read.table(x,header = T)})
  
  
  # create a GRangesList
  grs = lapply(segs, function(x) {GRanges(seqnames = Rle(x$chromosome),
                                          ranges = IRanges(start = x$start.pos,
                                                           end = x$end.pos),
                                          mcols = data.frame(copy_num = as.character(x$CNt),
                                                             A = as.character(x$A),
                                                             B = as.character(x$B)))})
  names(grs) = paste0("tumor",1:length(grs))
  grl = GRangesList(grs)
  
  
  intersected = Reduce(intersect, grs)
  intersected = intersected %>% 
    as.data.frame() %>% 
    split(.,.$seqname)
  
  if (sum(unlist(lapply(intersected, function(x) {nrow(x)})))==0) {
    print("There's no intersected copy neutral regions!")
   out = tibble(chrom = character(),
                chromStart = numeric(),
                chromEnd = numeric()) 
  } else {
    out = Reduce(bind_rows,intersected)
    
    out = out %>% 
      rename(chrom = seqnames,
             chromStart = start,
             chromEnd = end) %>% 
      select(chrom,chromStart,chromEnd)
  }
  
  num_copyneutral_segs = bind_rows(num_copyneutral_segs,
                                   tibble(Unified_ID = pat,
                                          num_copyneutral_segs = nrow(out)))
  
  write.table(out,paste0(odir,"/",pat,"_CopyNeutralSegments.bed"),
              sep = "\t",row.names = F,quote = F)
}

write_tsv(num_copyneutral_segs,paste0(odir,"/num_copyneutral_segs.tsv"))
