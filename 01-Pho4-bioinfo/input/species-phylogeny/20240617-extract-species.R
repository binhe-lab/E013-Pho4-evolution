# title: extract species for plotting a species tree
# author: Bin He
# date: 2023-07-07, modified 2024-06-17

# Load libraries
require(tidyverse)
require(ggtree)
require(treeio)

# Use TreeHouse to extract a subset of species
## running the app
shiny::runGitHub("treehouse", "JLSteenwyk")
## save and edit outside this script

# Read tree associated info
spsInfo <- read_tsv("fig1-species-phylogeny/20240619-subset-20-species-info.tsv", col_types = cols())

# Read ML tree with bootstrap values
sps.tree.ML <- read.tree("fig1-species-phylogeny/20240618-subset-20-ML-bootstrap-species-tree-reroot") %>% 
  as_tibble() %>% 
  #mutate(label = gsub("_", " ", label)) %>% 
  left_join(spsInfo, by = c("label" = "treeName")) %>% 
  as.treedata() #%>% 

# Read time calibrated tree
sps.tree.time <- read.tree("fig1-species-phylogeny/20240618-subset-20-time-calibrated-species-tree") %>% 
  as_tibble() %>% 
  #mutate(label = gsub("_", " ", label)) %>% 
  left_join(spsInfo, by = c("label" = "treeName")) %>% 
  as.treedata() #%>% 
  #root(outgroup = "Yarrowia_lipolytica", edgelabel = TRUE)

# Plot the tree
p.tree.ML <- ggtree(sps.tree.ML, ladderize = TRUE) + #scale_y_reverse() + #xlim(0,3) +
  geom_tiplab(aes(label = species), size = 10, face = 3, as_ylab = TRUE) +
  #geom_tiplab(size = 3.2, fontface = "italic", align = TRUE, linesize = 0.1, offset = 0.05) +
  geom_treescale(x = 0.1, width = 0.5, linesize = 1.2) +
  #geom_hilight(node = clade["MDR"], fill = "#7F00FF", alpha = 0.15)  + # MDR
  #geom_hilight(node = clade["CaLo"], fill = "pink", alpha = 0.25)    + # Candida/Lodderomyces
  #geom_hilight(node = clade["glabrata"], fill = "steelblue", alpha = 0.15)  + # glabrata
  geom_tippoint(aes(color = pathogen), size = 2) +
  scale_color_manual(values =  c("crustacean" = "#6a5acd",
                                 "human" = "#d14949", 
                                 "human (rare)" = "steelblue",
                                 "no report" = "gray20")) +
  #guides(color = guide_legend(byrow = TRUE)) +
  theme(legend.position = c(0.227, 0.90))
ggsave(p.tree.ML, file = "fig1-species-phylogeny/20240619-species-tree-ML-bootstrap.png", width = 4, height = 5)

p.tree.time <- ggtree(sps.tree.time, ladderize = TRUE) + #scale_y_reverse() + #xlim(0,3) +
  theme_tree2() +
  geom_tiplab(aes(label = species), size = 10, face = 3, as_ylab = TRUE) +
  #geom_tiplab(size = 3.2, fontface = "italic", align = TRUE, linesize = 0.1, offset = 0.05) +
  #geom_treescale(x = 0.1, width = 1, linesize = 1.2) +
  #geom_hilight(node = clade["MDR"], fill = "#7F00FF", alpha = 0.15)  + # MDR
  #geom_hilight(node = clade["CaLo"], fill = "pink", alpha = 0.25)    + # Candida/Lodderomyces
  #geom_hilight(node = clade["glabrata"], fill = "steelblue", alpha = 0.15)  + # glabrata
  geom_tippoint(aes(color = pathogen), size = 2) +
  scale_color_manual(values =  c("crustacean" = "#6a5acd",
                                 "human" = "#d14949", 
                                 "human (rare)" = "steelblue",
                                 "no report" = "gray20")) +
  #guides(color = guide_legend(byrow = TRUE)) +
  theme(legend.position = c(0.227, 0.90))
ggsave(p.tree.time, file = "fig1-species-phylogeny/20240619-species-tree-time-calibrated.pdf", width = 4, height = 5)
