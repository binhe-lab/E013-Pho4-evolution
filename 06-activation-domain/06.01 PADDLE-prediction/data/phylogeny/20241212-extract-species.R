# title: extract species for plotting a species tree
# author: Bin He
# date: 2023-07-07, modified 2024-12-12

# Load libraries
require(tidyverse)
require(ggtree)
require(treeio)

# Use TreeHouse to extract a subset of species
## running the app
## ---
## uncomment the following line to run the script
## shiny::runGitHub("treehouse", "JLSteenwyk")
## ---
## save and edit outside this script

# Read tree associated info
spsInfo <- read_tsv("20241212-sub-8-species-info.tsv", col_types = cols())

# Determine the species order
#sps.order <- c("C. albicans", "L. kluyveri", "N. castellii", "N. bacillisporus",
#               "C. bracarensis", "C. glabrata", "S. mikatae", "S. cerevisiae")

# Read time calibrated tree
sps.tree.time <- read.tree("20241212-sub-8-species-time-calibrated-tree.nwk") %>% 
  #rotateConstr(sps.order) %>% 
  as_tibble() %>% 
  left_join(spsInfo, by = c("label" = "treeName")) %>% 
  as.treedata()

p.tree.time <- sps.tree.time %>% 
  ggtree(ladderize = FALSE) + #scale_y_reverse() + #xlim(0,3) +
  theme_tree2() +
  geom_tiplab(aes(label = species), size = 16, face = 3, as_ylab = TRUE) +
  geom_tippoint(size = 2)
  #geom_tiplab(size = 3.2, fontface = "italic", align = TRUE, linesize = 0.1, offset = 0.05) +
  #geom_treescale(x = 0.1, width = 1, linesize = 1.2) +
  #geom_hilight(node = clade["MDR"], fill = "#7F00FF", alpha = 0.15)  + # MDR
  #geom_hilight(node = clade["CaLo"], fill = "pink", alpha = 0.25)    + # Candida/Lodderomyces
  #geom_hilight(node = clade["glabrata"], fill = "steelblue", alpha = 0.15)  + # glabrata
ggsave(p.tree.time, file = "20241212-species-tree-time-calibrated.png", width = 3, height = 5)
