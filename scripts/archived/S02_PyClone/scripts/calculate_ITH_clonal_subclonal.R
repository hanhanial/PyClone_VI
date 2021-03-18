library(tidyverse)

# This function removes singleton cluster and calculate the number of subclonal mutations to 
# the number of clonal mutations.

args = commandArgs(trailingOnly=TRUE)

patient = args[1]

# if (file.exists("./SC_to_C_ratio.txt")) {
#   file.remove("./SC_to_C_ratio.txt")
# }

# Function to do backmapping
rearrange_int <- function(a){
  b <- c()
  for (i in seq(1, length(a))){
    if (a[i] == i){
      b[i] = i
    } else{
      b[a[i]] = i
    }
  }
  return(b)
}

loci = read_tsv(paste0("./tables/loci.tsv"), col_types = cols())
cluster = read_tsv(paste0("./tables/cluster.tsv"), col_types = cols())

# Filter out singleton clusters
cluster = cluster %>% filter(! size == 1)
loci = loci %>% filter(cluster_id %in% cluster$cluster_id)

# Start cluster ID at 1 since stupid R uses 1-indexing...
cluster$cluster_id <- cluster$cluster_id + 1
loci$cluster_id <- loci$cluster_id + 1

# Rearrange cluster id to make it contiguous for easier downstream analysis
clust_mapID <- rearrange_int(unique(cluster$cluster_id))
cluster$cluster_id <- clust_mapID[cluster$cluster_id]
loci$cluster_id <- clust_mapID[loci$cluster_id]

# Write filtered loci and cluster TSV
loci = write_tsv(loci, paste0(patient, "_loci_filtered.tsv"))
cluster = write_tsv(cluster, paste0(patient, "_cluster_filtered.tsv"))

spread_clust <- spread(cluster %>% select(-std, -size), "sample_id", "mean")

write_tsv(spread_clust %>% select(-cluster_id), paste0(patient, "_cluster_prev.txt"), col_names = F)
write_tsv(as.data.frame(colnames(spread_clust %>% select(-cluster_id))),
	  paste0(patient, "_sample_order.txt"),
	  col_names=F)

clust_averageCCF = cluster %>% group_by(cluster_id) %>% summarise(meanCCF = mean(mean))
clust_maxCCF = clust_averageCCF[which.max(clust_averageCCF$meanCCF), 'cluster_id', drop=T]

# Visualize CCF across samples for cluster
ggplot(data=cluster) + geom_line(aes(x=sample_id, y=mean, group=cluster_id, color=factor(cluster_id)))
ggsave(paste0(patient, "_clone_prevalence.pdf"), width=16, height=12)

# Get proportion of mutations in highest CCF cluster vs the rest
numClonal = loci %>% filter(cluster_id == clust_maxCCF) %>% select(mutation_id) %>% unique %>%  nrow()
numSubclonal = loci %>% filter(!cluster_id == clust_maxCCF) %>% select(mutation_id) %>% unique %>%  nrow()
SC_to_C_ratio = numSubclonal/numClonal
write(paste0("Patient ", patient, " subclonal to clonal driver ratio: ", SC_to_C_ratio),
      "SC_to_C_ratio.txt", append = T)
