# ---
# title: functions for Pho4 chimera data plotter v3
# author: Bin He
# date: 2023-11-29
# ---

# Define functions
# 1. Data selection
# given a pattern, select chimera by symbols matching the pattern
my_data_select <- function(pattern = NULL, Set = NULL){
  # function to select a subset of the chimeras given a set of rules
  # validate the input
  
  # change region 4 into a consistent format
  tmp <- ximera %>% 
    mutate(
      symbol = as.character(symbol),
      Symbol = ifelse(
        nchar(symbol) == 5,
        paste0(str_sub(symbol, 1, 3), 
               str_sub(symbol, 4, 4), 
               str_sub(symbol, 4, 4), 
               str_sub(symbol, 5, 5)),
        symbol
      )) %>% 
    select(plasmid, Symbol)
  # starting set
  if(length(Set) == 0)
    xim <- filter(tmp, !plasmid %in% refs)
  else
    xim <- filter(tmp, set %in% Set, !plasmid %in% refs)
  # compare to the pattern
  try(if(nchar(pattern) != 6) stop("Pattern must be a string with 6 characters"))
  symbols = xim$Symbol # extrac the symbols for testing
  include = nchar(symbols) > 0 # initialize the inclusion vector
  for(i in 1:6){
    p = substr(pattern, i, i)
    if(p != "X" & p != "x"){ # ignore X and x
      test = toupper(str_sub(symbols, i, i)) == toupper(p)
      include = include & test
    }
  }
  select <- cbind(xim, include)
  return(select$plasmid[select$include])
}

# combine sets
my_data_select_m <- function(patterns = NULL, Set = NULL){
  # combine multiple selections
  all_selected = c()
  try(if(length(patterns) == 0) stop("No patterns provided"))
  for(i in patterns)
    all_selected = c(all_selected, my_data_select(pattern = i))
  return(unique(all_selected))
}

# not in use
my_data_add <- function(input){
  # given a comma separated set of symbols, return the corresponding plasmid IDs
  symbols <- unlist(str_split(input, pattern = "\\s*,\\s*"))
  ximera$plasmid[ximera$symbol %in% symbols]
}

# combine two sets
my_selection_combine <- function(set1, set2){
  # given two vectors, combine them and return the unique set of plasmid IDs
  unique(c(set1, set2))
}

# 2. Data prep and transform
my_data_prep <- function(selection){
  # given a selection of chimera ID (plasmid), prepare a data frame for plotting
  # subset data
  tmp <- ximera %>% 
    filter(plasmid %in% c(refs, selection)) %>% 
    select(plasmid, symbol, group) %>% 
    inner_join(dat, by = "plasmid") %>% 
  return(tmp)
}

# 3. Plotting
# common theme elements
themes <- list(
  theme_cowplot(font_size = 20),
  panel_border(color = "gray30"),
  background_grid(major = "y", minor = "none"),
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "courier"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_markdown())
)

# plot R/G ratios and the individual components
my_plot_ratio <- function(selection){
  # plot R/G ratios and the individual components
  # this is the "raw data plot" for examining the underlying evidence
  # for higher level plots like the relative activity below
  tmp <- my_data_prep(selection)
  p <- tmp %>% 
    select(-c(FSC.H, nGFP, nRFP, flag)) %>% 
    mutate(`R/G` = YL2.H/BL1.H) %>% 
    pivot_longer(cols = c(BL1.H, YL2.H, `R/G`), 
                 names_to = "parameter", values_to = "value") %>% 
    mutate(parameter = factor(parameter, levels = c("R/G", "YL2.H", "BL1.H"),
                              labels = c("RFP/GFP", "PHO5pRFP", "Pho4-GFP"))) %>% 
    ggplot(aes(x = symbol, y = value, group = host)) + 
    stat_summary(aes(group = host), fun.data = "mean_cl_boot", geom = "errorbar",
                 position = position_dodge(0.6), width = 0.3) +
    geom_bar(aes(fill = host), alpha = 0.8, width = 0.6,
             stat = "summary", fun = "mean", position = position_dodge(0.6)) +
    geom_point(data = function(x) subset(x, !symbol %in% c("CCCCC", "SSSSS")),
               aes(group = host, color = host), size = 1, shape = 3, alpha = 0.9,
               position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1)) +
    scale_color_manual(values = point.colors) +
    scale_fill_manual(values = host.colors) +
    facet_grid(parameter~group, scales = "free", space = "free_x", switch = "y") +
    xlab("Pho4 chimera") + themes
  return(p) 
}

my_plot_boost <- function(selection){
  # given a selection of chimera IDs, plot their functionality w/PHO2
  # relative to ScPho4, and their boost
  # dat
  tmp <- filter(ximera, plasmid %in% c(refs, selection)) %>% 
    mutate(perc_pho2 = A_pho2/A_PHO2) %>% 
    pivot_longer(cols = c(s_PHO2, boost, perc_pho2), 
                 names_to = "parameter", values_to = "ratio")
  # labeller
  par.explain <- c(
    s_PHO2 = "Rel. A<sub>PHO2</sub>",
    boost = "Boost",
    perc_pho2 = "%A<sub>pho2∆</sub>"
  )
  p <- ggplot(tmp, aes(x = symbol, y = ratio)) +
    geom_col(width = 0.5, color = "black", fill = "gray80") +
    geom_hline(yintercept = 1, linetype = 2, color = "gray30") +
    facet_grid(parameter~group, scales = "free", space = "free_x", switch = "y",
              labeller = labeller(parameter = par.explain)) +
    #theme_bw(base_size = 18) +
    #background_grid(minor = "none") +
    themes
  return(p)
}

my_scatter_plot <- function(pattern){
  # plot all chimeras by their A_PHO2 and A_pho2∆, highlighting a subset of
  # the chimeras by their composition.
  selection = my_data_select(pattern = pattern)
  scatter.colors = c("ScPho4" = "forestgreen", "CgPho4" = "blue3", "cyan2", "other" = "gray20")
  names(scatter.colors)[3] = pattern
  p <- ximera %>% 
    mutate(A_PHO2 = signif(A_PHO2, digits = 2),
           A_pho2 = signif(A_pho2, digits = 2),
           group = case_when(
             symbol == "CCCCC" ~ "CgPho4",
             symbol == "SSSSS" ~ "ScPho4",
             plasmid %in% selection ~ pattern,
             .default = "other"
           ),
           group = fct_relevel(group, names(scatter.colors))) %>% 
    ggplot(aes(x = A_PHO2, y = A_pho2, label = symbol)) + 
    geom_point(aes(color = group), size = 2.5) + 
    scale_color_manual("Pho4 type", values = scatter.colors) +
    geom_abline(slope = 1) +
    theme_gray(base_size = 14)
  return(p)
}

