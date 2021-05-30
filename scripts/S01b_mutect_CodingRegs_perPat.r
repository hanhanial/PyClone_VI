'
Script: S01b
for each patient, combine coding regions from their samples 
--> those merged regions will be used as bed file to filter mutations
'

suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))


merged_segs = function(segs) { # segs is the a vector of directories to segments 
  # list of segment data of all the tumor samples of the patient
  dat = lapply(segs, function(x) {dat = read.table(x,header = T)})
  
  
  # create a GRangesList
  grs = lapply(dat, function(x) {GRanges(seqnames = Rle(x$chrom),
                                         ranges = IRanges(start = x$chromStart,
                                                          end = x$chromEnd))})
  names(grs) = paste0("tumor",1:length(grs))
  grl = GRangesList(grs)
  
  
  merged = Reduce(union, grs)
  merged = merged %>% 
    as.data.frame() %>% 
    split(.,.$seqname)
  
  out = Reduce(bind_rows,merged)
  out = out %>% 
    rename(chrom = seqnames,
           chromStart = start,
           chromEnd = end) %>% 
    filter(!grepl(pattern = "_",x = chrom)) %>% 
    select(chrom,chromStart,chromEnd)
  
  return(out)
}



args = commandArgs(trailingOnly = T)

# master file, or list of DNA libs to be processed for each patient
mf = args[1]

# working directory
wdir = args[2]

'
mf = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/pyclone_masterfile.csv"
wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210429/"
'

mf = read_csv(mf)

# directory storing bed files (SNP and coding muts) of each sample
idir = paste0(wdir,"/S01a_mutect_CodingRegs_perSam/")

# output directory to store bed files (SNP and coding muts) of each patient
odir = paste0(wdir,"/S01b_mutect_CodingRegs_perPat/")
system(paste0("mkdir -p ",odir))


# DNA libs and directory to the file that storing their copy-neutral segments
snp_segs = list.files(idir,pattern = "_SNPs.bed")
snp_segs = tibble(DNA_lib = snp_segs,
                  snp_seg = paste0(idir,snp_segs)) %>% 
  mutate(DNA_lib = gsub(pattern = "_SNPs.bed",replacement = "",x = DNA_lib))

coding_segs = list.files(idir,pattern = "_coding.bed")
coding_segs = tibble(DNA_lib = coding_segs,
                     coding_seg = paste0(idir,coding_segs)) %>% 
  mutate(DNA_lib = gsub(pattern = "_coding.bed",replacement = "",x = DNA_lib))

all_segs = inner_join(snp_segs,coding_segs)
all_segs = all_segs %>% 
  filter(DNA_lib %in% mf$DNA_lib) %>% 
  inner_join(mf)

pats = mf %>% 
  pull(Unified_ID) %>% 
  unique()

# to record number of (merged) SNP/coding SNP for each patient
num_segs = tibble(Unified_ID = character(),
                  num_SNP = numeric(),
                  num_coding = numeric())

for (p in 1:length(pats))
{
  # p = 1
  pat = pats[p]
  
  print(paste0("Processing patient ",p,": ",pat))
  
  pat_snps = all_segs %>% 
    filter(Unified_ID==pat) %>%
    pull(snp_seg)
  pat_merged_snps = merged_segs(segs = pat_snps)
  
  pat_coding = all_segs %>% 
    filter(Unified_ID==pat) %>%
    pull(coding_seg)
  pat_merged_coding = merged_segs(segs = pat_coding)
  
  num_segs = bind_rows(num_segs,
                       tibble(Unified_ID = pat,
                              num_SNP = nrow(pat_merged_snps),
                              num_coding = nrow(pat_merged_coding)))
  
  write.table(pat_merged_snps,
              paste0(odir,"/",pat,"_SNPSegments.bed"),
              sep = "\t",row.names = F,quote = F)
  
  write.table(pat_merged_coding,
              paste0(odir,"/",pat,"_CodingSegments.bed"),
              sep = "\t",row.names = F,quote = F)
}

write_tsv(num_segs,paste0(odir,"/num_merged_segs.tsv"))


