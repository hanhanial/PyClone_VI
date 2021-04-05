'
script to get 2 sets of unique mutations (coding or non-coding) for each patient
--> from the unique sets, randomly select as many coding as we can, then top up with non-coding mutations
--> merge with purity info
--> pyclone input for each patient
'

suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly = T)

# patient ID
pat = args[1]

# number of muts to be input to pyclone
num_muts = args[2]
num_muts = as.numeric(num_muts)

# random seed for selection mutations
random_seed = args[3]

# working dir
wdir = args[4]

# input dir
idir = paste0(wdir,"/assign_CN/",pat,"_pyclone_input/")


'
pat = "E001"
num_muts = 3000
random_seed = 123
wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210318//S04_pyclone"
idir = paste0(wdir,"/assign_CN/",pat,"_pyclone_input/")

'

# output dir 
odir = pat
system(paste0("mkdir -p ",odir))

mut_files = list.files(idir)

# patient's coding files
coding_files = mut_files[startsWith(x = mut_files,prefix = "coding")]
coding_files = tibble(sample_id = gsub(pattern = "coding_",replacement = "",x = coding_files),
                      file = paste0(idir,coding_files)) %>% 
  mutate(sample_id = gsub(pattern = ".tsv",replacement = "",x = sample_id))
coding_files = lapply(1:nrow(coding_files), function(x) {
  dat = read_tsv(coding_files$file[x]) %>% 
    mutate(sample_id = as.character(coding_files$sample_id[x]))
  return(dat)
})

# patient's non-coding files
noncoding_files = mut_files[startsWith(x = mut_files,prefix = "noncoding")]
noncoding_files = tibble(sample_id = gsub(pattern = "noncoding_",replacement = "",x = noncoding_files),
                      file = paste0(idir,noncoding_files)) %>% 
  mutate(sample_id = gsub(pattern = ".tsv",replacement = "",x = sample_id))
noncoding_files = lapply(1:nrow(noncoding_files), function(x) {
  dat = read_tsv(noncoding_files$file[x]) %>% 
    mutate(sample_id = as.character(noncoding_files$sample_id[x]))
  return(dat)
})



# list of unique coding muts
coding_muts = lapply(coding_files, function(x) {
  return(x$mutation_id)
})
coding_muts = unique(Reduce(f = c,x = coding_muts))


# list of unique non-coding muts
noncoding_muts = lapply(noncoding_files, function(x) {
  return(x$mutation_id)
})
noncoding_muts = unique(Reduce(f = c,x = noncoding_muts))


# randomly pick coding muts 
# - if length(coding_muts)<=num_muts, use all coding muts
# - otherwise (i.e. more coding muts than num muts to be selected), randomly select num_muts from coding_muts
if (num_muts>=length(coding_muts)) {
  selected_coding_muts = coding_muts
} else { 
  set.seed(random_seed)
  selected_coding_muts = sample(x = coding_muts,size = num_muts,replace = F) }


# get the list of all selected muts 
# - if we have enough coding muts (i.e. length(selected_coding_muts)>=num_muts), then all selected muts will be coding muts
# - otherwise, use noncoding_muts to top up (num_muts-length(coding_muts)) muts 
if (length(selected_coding_muts)>=num_muts) {
  selected_muts = selected_coding_muts
} else {
  num_selected_noncoding_muts = num_muts-length(coding_muts)
  
  if (num_selected_noncoding_muts>=length(noncoding_muts)) {
    # if we dont have enough non-coding muts, use all of them
    selected_noncoding_muts = noncoding_muts
  } else {
    # else randomly select num_selected_noncoding_muts muts from noncoding_muts
    set.seed(random_seed)
    selected_noncoding_muts = sample(x = noncoding_muts,size = num_selected_noncoding_muts,replace = F)
  }
  
  selected_muts = c(selected_coding_muts,selected_noncoding_muts)
}


# pyclone input
coding_inp = lapply(coding_files, function(x) {
  x %>% filter(mutation_id %in% selected_muts)
})
coding_inp = Reduce(f = bind_rows,x = coding_inp)

noncoding_inp = lapply(noncoding_files, function(x) {
  x %>% filter(mutation_id %in% selected_muts)
})
noncoding_inp = Reduce(f = bind_rows,x = noncoding_inp)

pat_inp = bind_rows(coding_inp,noncoding_inp) %>% 
  arrange(sample_id)

# duplicated muts
duplicated_muts = pat_inp[duplicated(pat_inp[,c("sample_id","mutation_id")]),] %>% 
  select(sample_id,mutation_id,everything()) %>% 
  arrange(sample_id,mutation_id)
if (nrow(duplicated_muts)>0) {
  write_tsv(duplicated_muts,paste0(pat,"/duplicated-mutations.tsv"))
}


# Purity file
purity = paste0(wdir,"/seg_purity/",pat,".purity")
purity = read_tsv(purity)
purity = purity %>% 
  rename(sample_id = Tumor_Sample_Barcode,
         tumour_content = purity) %>% 
  select(sample_id,tumour_content)


# merge to get tumor content
pat_inp = left_join(pat_inp,purity) 

pat_inp1 = pat_inp %>% 
  rename(alt_counts = var_counts) %>% 
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
write_tsv(pat_inp1,paste0(pat,"/",pat,".tsv"))

# output number of selected coding/noncoding muts
pat_mut_count = tibble(Unified_ID = pat,
                       num_muts = length(selected_muts),
                       num_coding = sum(coding_muts %in% selected_muts),
                       num_noncoding = sum(noncoding_muts %in% selected_muts))
write_tsv(pat_mut_count,paste0(pat,"/num_muts.tsv"))
                                              