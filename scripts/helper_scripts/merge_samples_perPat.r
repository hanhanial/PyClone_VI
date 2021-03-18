suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly = T)

# patient ID
# pat_id = "PB_16_1537"
pat_id = args[1]

# working dir
# wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/output/S01_prepareInputs/S05_pyclone_inputs"
wdir = args[2]

# Purity file
purity = paste0(wdir,"/seg_purity/",pat_id,".purity")
purity = read_tsv(purity)
purity = purity %>% 
  rename(sample_id = Tumor_Sample_Barcode,
         tumour_content = purity) %>% 
  select(sample_id,tumour_content)


# list of patient's samples
idir = paste0(wdir,"/assign_CN/",pat_id,"_pyclone_input/")
pat_files = list.files(idir)
dat = lapply(pat_files, function(x) {
  # x = pat_files[1]
  o = read_tsv(paste0(idir,x))
  dna_lib = gsub(pattern = ".tsv",replacement = "",x)
  o = o %>% 
    rename(alt_counts =  var_counts) %>% 
    mutate(sample_id = dna_lib)
})

dat1 = Reduce(f = bind_rows,x = dat)
# dat1[duplicated(dat1),]

# get tumor content
dat1 = left_join(dat1,purity)

dat2 = dat1 %>% 
  select(mutation_id,sample_id,
         ref_counts,alt_counts,
         normal_cn,major_cn,minor_cn,
         tumour_content) %>% 
  group_by(mutation_id,sample_id) %>% 
  summarise_all(.funs = mean) %>% 
  mutate_at(vars(ref_counts,alt_counts,
                 normal_cn,major_cn,minor_cn), 
            funs(round)) %>% 
  ungroup()


write_tsv(dat2,paste0(pat_id,".tsv"))