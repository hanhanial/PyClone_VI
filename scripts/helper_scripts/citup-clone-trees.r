## A set of common functions I use across other scripts

# Colormap from stackexchange, good for class palette
c25 = c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

# Coding definition (non-syn including silent) for TCGA MAF format
coding_def = c("De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
                "Start_Codon_Del", "Start_Codon_SNP", "Frame_Shift_Del",
                "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                "Silent", "Splice_Site", "Start_Codon_Ins",
                "Translation_Start_Site")


'
- patient = patient ID
- citup_dir = directory to citup result of the patient
- mf = master file, which maps DNA_lib to Tumor
- nPoints = number of segments used to draw clonal decomposition for each tumor sample
higher nPoints --> more defined decomposition, but also output much bigger size file/pdf
'
clone_tree = function(patient,citup_dir,mf,nPoints) {
  
  #### Clone tree of the patient ####
  # Read in the assignment of clusters of mutations to different clone IDs by citup
  var_assignment = read.delim2(list.files(path=citup_dir, 
                                          pattern=".*variant_assignment.*.txt$", 
                                          recursive = FALSE, full.names = TRUE),
                               header = FALSE)
  var_assignment = var_assignment %>% 
    rename(clone_id = V1) %>% 
    mutate(cluster_id = 1:nrow(var_assignment))
  
  # Read in list of mutations for each cluster
  mutation = read.delim2(list.files(path=citup_dir, 
                                    pattern=".*cluster_filtered.*.tsv$", 
                                    recursive = FALSE, full.names = TRUE))
  
  mutation = mutation %>% 
    left_join(var_assignment)
  
  # get mean CCF of each clone
  clone_CCF = mutation %>% 
    mutate(cellular_prevalence = as.numeric(as.character(cellular_prevalence))) %>% 
    group_by(clone_id) %>% 
    summarise(CCF = mean(cellular_prevalence,na.rm = T))
  
  # get number of mutations of each clone
  clone_size = mutation %>% 
    distinct(clone_id,cluster_id,cluster_size) %>% 
    group_by(clone_id) %>% 
    summarise(num_muts = sum(cluster_size))
  
  
  
  # Read in adjacency list (i.e. which clone is linked to which clone on the tree?)
  adjacency = read.delim2(list.files(path=citup_dir, 
                                     pattern=".*tree.*.adj$", 
                                     recursive = FALSE, full.names = TRUE),
                          header = FALSE)
  
  
  # For igraph, edge starts from 1. Citup output starts from 0.
  # I manually add 2 so it starts from 2. Then add a pseudo root node 1 to draw a
  # line to indicate root.
  adjacency = rbind(c(-1,0), adjacency)
  adjacency_igraph = graph_from_edgelist(as.matrix(adjacency) + 2)
  
  # Set vertices' names
  V(adjacency_igraph)$Clones = as.character(min(adjacency):max(adjacency))
  
  # Set root node name to none
  V(adjacency_igraph)$Clones[1] = ""
  
  # number of mutations in each clone
  V(adjacency_igraph)$num_muts = clone_size$num_muts[match(V(adjacency_igraph)$Clones, 
                                                           clone_size$clone_id)]
  V(adjacency_igraph)$num_muts[is.na(V(adjacency_igraph)$num_muts)] = 0
  
  # mean CCF in each clone
  V(adjacency_igraph)$CCF = clone_CCF$CCF[match(V(adjacency_igraph)$Clones, 
                                                clone_CCF$clone_id)]
  V(adjacency_igraph)$CCF[is.na(V(adjacency_igraph)$CCF)] = 0
  
  adjacency_tbl = as_tbl_graph(adjacency_igraph)
  
  graph_plot = ggraph(adjacency_tbl, 'tree') +
    
    # Use white for the pseudo root node to hide it.
    scale_colour_manual(values = c("#FFFFFF", c25[2:length(V(adjacency_igraph))])) +
    geom_edge_link(hjust=0.5, vjust=1.2, show.legend = FALSE,edge_width = 1.5,
                   label_alpha=1,label_size = 2.5, angle_calc = "along") +
    
    # size of the nodes (clones) is proportional to number of muts in that clone
    geom_node_point(aes(size=num_muts, 
                        color=Clones), 
                    show.legend = c(size=FALSE,
                                    color=TRUE,
                                    alpha=FALSE)) +
    # Scale point to make bigger as some nodes don't have a lot of muts
    scale_size_continuous(range = c(5,15)) +
    
    ggtitle(patient) +
    theme_graph(foreground = 'steelblue', fg_text_colour = 'white',
                base_family = "sans") +
    theme(plot.title = element_text(size = 18,hjust = 0.5),
          panel.border = element_blank()) +
    scale_y_reverse() 
  
  # Get x and y range to give space to labels
  xlims = ggplot_build(graph_plot)$layout$panel_scales_x[[1]]$range$range
  xlims[1] = xlims[1] - 0.25
  xlims[2] = xlims[2] + 0.25
  
  graph_plot = graph_plot + xlim(xlims)
  
  
  
  
  #### Clonal decomposition for each tumor sample ####
  # get sample order
  sample_order = read_tsv(list.files(path=citup_dir, 
                                     pattern=".*sample_order*.txt$", 
                                     recursive = FALSE, full.names = TRUE),
                          col_names = FALSE)
  # prevalence of each clone in each sample
  clone_freq = read_tsv(list.files(path=citup_dir, 
                                   pattern=".*clone_freq*.txt$", 
                                   recursive = FALSE, full.names = TRUE))
  
  clone_freq$sample = sample_order[clone_freq$sample + 1, "X1", drop=TRUE]
  
  plot = NULL
  plot[["graph"]] = graph_plot
  for (sector in unique(clone_freq$sample))
  {
    tumor = mf$Tumor[which(mf$DNA_lib==sector)]
    
    # Get radial position for labelling
    sector_freq = clone_freq %>% 
      filter(sample == !!sector) %>% 
      arrange(clone_id) %>%
      mutate(rad_pos = cumsum(clonal_prev) - clonal_prev/2) %>% ungroup()
    
    sector_freq = sector_freq %>% 
      mutate(clone_label = if_else(clonal_prev > 0.01, 
                                   as.character(round(clonal_prev, 2)), ""))
    
    color_clone = c25[sector_freq$clone_id + 2]
    names(color_clone) = sector_freq$clone_id
    sector_freq$clone_id = factor(sector_freq$clone_id, levels=rev(sector_freq$clone_id))
    
    # Generate random circle
    angles = runif(nPoints) * 2 * pi
    rads = 1 * sqrt(runif(nPoints))
    
    # Generate cartesian coord
    x = rads * cos(angles)
    y = rads * sin(angles)
    
    df = tibble(x=x, y=y, distance = rads)
    
    
    # Get number of points for ccf
    sector_freq = sector_freq %>%
      mutate(clonal_prev_integer = floor(clonal_prev * nPoints))
    
    # Calculate distance matrix for the "tumor" points
    dist_mat = as.matrix(dist(df))
    
    # Run through each clone and select points closest
    # to a randomly sampled centroid
    rand_cen = NULL
    selected = NULL
    
    # remove clones with too few points
    sector_freq = sector_freq[which(sector_freq$clonal_prev_integer >= floor(0.01 * nPoints)),] #floor(0.01 * nPoints)
    
    # Start growing at the edge
    rand_i = which.max(df$distance)
    
    for (i in 1:nrow(sector_freq))
    {
      # Change to character to index by rownames
      rand_i = as.character(rand_i)
      
      # Select number of X points closest to the rand selected index
      sorted_dist = names(sort(dist_mat[rand_i, ]))
      closest_points = sorted_dist[1:(sector_freq[i, ]$clonal_prev_integer)]
      
      # Select next closest point as the next clone to grow
      rand_i = sorted_dist[(sector_freq[i, ]$clonal_prev_integer + 1)]
      
      # Select the closest point + the randomly sampled index
      selected = unique(c(selected, closest_points))
      
      # Removed selected points from distance matrix
      dist_mat = dist_mat[!rownames(dist_mat) %in% selected, !colnames(dist_mat) %in% selected]
      
      df[closest_points, "clone_id"] = sector_freq[i, "clone_id"]
    }
    
    if (sum(is.na(df$clone_id)) == 0) {
      df1 = df
    } else {
      df1 = df[-which(is.na(df$clone_id)),]
    }
    
    plot[[tumor]] = ggplot(df1, aes(x = x, y = y)) +
      
      geom_voronoi_tile(aes(x = x,y = y,fill = clone_id, group=1),alpha=0.8) +
      scale_fill_manual(values = color_clone) + 
      
      geom_voronoi_segment(aes(x = x,y = y, group=1),alpha=0.2) +
      
      xlim(c(-1.0, 1.0)) +
      
      ylim(c(-1.0, 1.0)) +
      coord_fixed() + 
      theme_void() + 
      theme(legend.position="none", plot.title = element_text(family = "sans", size = 20,hjust = 0.5)) +
      ggtitle(tumor)
    
    
  }
  
  return(plot)
}