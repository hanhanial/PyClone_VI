library(tidyverse)
library(ggraph)
library(tidygraph)
library(igraph)
library(gridExtra)
library(ggrepel)
library(ggpubr)
library(ggforce)

source("/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts/citup-clone-trees.r")


# master file
mf = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S21_PyClone_VI/d20210622/S01_input_masterfile.csv"

# working directory
wdir = "/mnt/projects/lailhh/workspace/Metastasis_Feb2020/S21_PyClone_VI/d20210622/S02_PyClone_VI/output/"


mf = read_csv(mf)
mf = mf %>% 
  select(Unified_ID,Tumor,DNA_lib)

# create output directory
odir = paste0(wdir,"/S05_citup-clone-trees/")
system(paste0("mkdir -p ",odir))


# Load drivers for filtering
drivers = load_drivers()


citup_dir = list.files(paste0(wdir,"/S04_pyclone/citup_results"))
citup_dir = tibble(Unified_ID = citup_dir) %>% 
  mutate(citup_dir = paste0(wdir,"/S04_pyclone/citup_results/",Unified_ID,"/"))


all_pat_plots = NULL
arranged_all_pat_plots = NULL
for (k in 1:nrow(citup_dir))
{
  # k = 1
  patient = citup_dir$Unified_ID[k]
  
  print(paste0("Processing patient ",k,": ",patient))
  
  pat_citup_dir = citup_dir$citup_dir[k]
  
  # Pyclone tree + clones for each tumor of the patient
  all_pat_plots[[patient]] = clone_tree(patient = patient,
                                        citup_dir = pat_citup_dir,
                                        mf = mf,nPoints = 180)
  
  
  # arrange all plots into a big plot
  plots = all_pat_plots[[patient]]
  numsec = length(plots) - 1
  layout_plot = rbind(c(rep(1, numsec)),
                      c(seq(2, numsec+1)))
  arranged_all_pat_plots[[patient]] = ggarrange(plots[[1]], 
                                                ggarrange(plotlist = plots[2:length(plots)], 
                                                          ncol=length(plots)-1, align="hv"), 
                                                nrow=2, heights = c(0.8,0.2))
}

saveRDS(all_pat_plots,paste0(odir,"all_pat_plots.rds"))
saveRDS(arranged_all_pat_plots,paste0(odir,"clonal_decomposition_plots.rds"))


pdf(paste0(odir,"clonal_decomposition_plots.pdf"),width = 24,height = 15)
for (j in 1:length(arranged_all_pat_plots))
{
  grid.draw(arranged_all_pat_plots[[j]])
}
dev.off()