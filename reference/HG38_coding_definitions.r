library(tidyverse)

vafs = readRDS("/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S08_trees/S02_DNAtrees/d20210223/S01_maf_per_Patient/VAF_perSam.rds")

variant_classifications = lapply(vafs, function(x) {
  return(unique(as.character(x$Variant_Classification)))
})

variant_classifications = unique(Reduce(f = c,x = variant_classifications))

(coding_defs = variant_classifications[!(variant_classifications %in% 
                                          c("RNA","IGR","Intron",
                                            "5'Flank","3'UTR",
                                            "5'UTR","COULD_NOT_DETERMINE"))])
(coding_defs = unique(c(coding_defs,"START_CODON_DEL")))

write_tsv(tibble(sort(coding_defs)),"HG38_coding_definitions.tsv",col_names = F)
