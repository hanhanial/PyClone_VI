'
- to prepare input for citup
- to visualize mean CCF of each cluster in each sample
'

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# pat_id = "A001"
pat_id = args[1]

# Function to do back-mapping
rearrange_int = function(a){
  b = c()
  for (i in seq(1, length(a))){
    if (a[i] == i){
      b[i] = i
    } else{
      b[a[i]] = i
    }
  }
  return(b)
}

# pyclone_res = read_tsv("/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/output/S01_prepareInputs/S05_pyclone_inputs/pyclone_results/A001/A001.tsv")
pyclone_res = read_tsv(paste0(pat_id,".tsv"))

# get correct names for sample_id and mutation_id
cluster = pyclone_res %>% 
  mutate(sample_id = gsub(pattern = "b'",replacement = "",x = sample_id),
         mutation_id = gsub(pattern = "b'",replacement = "",x = mutation_id)) %>% 
  mutate(sample_id = gsub(pattern = "'",replacement = "",x = sample_id),
         mutation_id = gsub(pattern = "'",replacement = "",x = mutation_id))

# get cluster size
cluster = cluster %>% 
  group_by(cluster_id) %>% 
  mutate(cluster_size = length(unique(mutation_id))) %>% 
  ungroup()

# remove cluster(s) with only one mutation in it
cluster = cluster %>% 
  filter(cluster_size!=1)

# Start cluster ID at 1 
cluster$cluster_id = cluster$cluster_id + 1

# Rearrange cluster id to make it contiguous for easier downstream analysis
clust_mapID = rearrange_int(unique(cluster$cluster_id))
cluster$cluster_id = clust_mapID[cluster$cluster_id]

# Write filtered cluster TSV
cluster = write_tsv(cluster, paste0(pat_id, "_cluster_filtered.tsv"))


cluster_meanCCF_perSam = cluster %>% 
  group_by(cluster_id,sample_id) %>% 
  summarise(meanCCF = mean(cellular_prevalence)) %>% 
  ungroup()
spread_clust = spread(data = cluster_meanCCF_perSam, 
                      key = "sample_id", value = "meanCCF")

write_tsv(x = spread_clust %>% select(-cluster_id), 
          paste0(pat_id,"_cluster_prev.txt"), 
          col_names = F)
write_tsv(x = as.data.frame(colnames(spread_clust %>% select(-cluster_id))),
          paste0(pat_id,"_sample_order.txt"),
          col_names=F)


# Visualize mean CCF across samples for cluster
if (length(unique(cluster_meanCCF_perSam$sample_id))==1) {
  ggplot(data=cluster_meanCCF_perSam) + 
    geom_point(aes(x=sample_id, y=meanCCF, 
                  group=cluster_id, color=factor(cluster_id)))
} else {
  ggplot(data=cluster_meanCCF_perSam) + 
    geom_line(aes(x=sample_id, y=meanCCF, 
                  group=cluster_id, color=factor(cluster_id)))
}

ggsave(paste0(pat_id, "_clone_prevalence.pdf"), width=16, height=12)
