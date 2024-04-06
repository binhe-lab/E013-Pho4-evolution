# ---
# title: functions for Pho4 chimera data plotter v3
# author: Bin He
# date: 2024-02-11
# ---

# Define plotting variables
# reference Pho4 plasmid ids
refs <- c("188", "194")
# colors
date.colors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))
point.colors = c("PHO2" = "gray10", "pho2" = "gray30")
host.colors = c("PHO2" = "gray50", "pho2" = "orange2")


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
    select(plasmid, Symbol, set)
  # starting set
  if(length(Set) == 0)
    xim <- filter(tmp, !plasmid %in% refs)
  else if (length(Set) == 1)
    xim <- filter(tmp, set == Set, !plasmid %in% refs)
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
    mutate(host = fct_recode(host, pho2 = "pho2∆"))
  return(tmp)
}

my_calc_region_effect <- function(region, variable){
  # this function takes the name of a variable of interest
  # x specifies the foreground region, which will be examined for its effect on
  # the variable of interest.
  # it then transforms the ximera data frame to preserve only the variable of
  # interest, pivots it wider after grouping by the background composition.
  
  # prepare the data by mutating the symbol column into fg and bg
  valid.var <- c("A_PHO2", "A_pho2", "rA_PHO2", "rA_pho2", "boost")
  if(!variable %in% valid.var)
    stop(paste0("Please specify one of the valid variable names:", 
                paste(valid.var, collapse = ", ")))
  tmp <- ximera %>% 
    filter(set == "M") %>% 
    select(plasmid, symbol, var = {{ variable }}) %>% 
    mutate(fg = str_sub(symbol, region, region) %>% toupper(),
           bg = symbol %>% toupper())
  # replace the foreground region with X for grouping
  str_sub(tmp$bg, region, region) <- "X"
  # reorganize the tibble for easier handling, optional
  tmp <- relocate(tmp, fg, bg, .before = symbol) %>% select(-symbol)
  # pivot the data into a wide format such that for each background, there
  # are two values for the variable of interest, one from the chimera with 
  # CgPho4's version in the foreground and another with ScPho4's version
  tmp <- tmp %>% 
    select(plasmid, fg, bg, var) %>% 
    pivot_wider(id_cols = bg, names_from = "fg", 
                values_from = c(plasmid, var)) %>% 
    unite(plasmid, starts_with("plasmid")) %>%
    mutate(label = paste(bg, plasmid, sep = "\n"))
  return(tmp)
}

my_comp_region_effect <- function(region){
  # this function uses my_calc_region_effect to get the value for the variable of interest
  # with either Cg or Sc version in the focal region, separately for each background composition
  # it does so for two variables, A_PHO2 and A_pho2, then calculate dA_PHO2, dA_pho2, and
  # combine them
  PHO2 = my_calc_region_effect(region, "A_PHO2") %>% 
    mutate(dA_PHO2 = var_C - var_S,
           # mean A_PHO2
           M_PHO2 = (var_S + var_C)/2,
           NF = ifelse(M_PHO2 <=3.5, TRUE, FALSE)) %>% 
    select(-var_S, -var_C)
  
  pho2 = my_calc_region_effect(region, "A_pho2") %>% 
    mutate(dA_pho2 = var_C - var_S, 
           M_pho2 = (var_S + var_C)/2) %>% 
    select(-var_S, -var_C)
  
  dat <- full_join(PHO2, pho2, by = c("bg", "plasmid", "label")) %>% 
    select(bg, plasmid, dA_PHO2, dA_pho2, M_PHO2, M_pho2, NF)
  
  return(dat)
}

# 3. Plotting
# common theme elements
themes <- list(
  theme_cowplot(),
  panel_border(color = "gray30"),
  background_grid(major = "y", minor = "none"),
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "courier"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_markdown())
)

# plot R/G ratios and the individual components
my_plot_components <- function(dat){
  # plot R/G ratios and the individual components
  # this is the "raw data plot" for examining the underlying evidence
  # for higher level plots like the relative activity below
  p <- dat %>% 
    select(-c(FSC.H, nGFP, nRFP, flag)) %>% 
    mutate(`R/G` = YL2.H/BL1.H) %>% 
    pivot_longer(cols = c(BL1.H, YL2.H, `R/G`), 
                 names_to = "parameter", values_to = "value") %>% 
    mutate(parameter = factor(parameter, levels = c("R/G", "YL2.H", "BL1.H"),
                              labels = c("RFP/GFP", "PHO5pRFP", "Pho4-GFP"))
    ) %>% 
    ggplot(aes(x = symbol, y = value, group = host)) + 
    stat_summary(aes(group = host), fun.data = "mean_cl_boot", geom = "errorbar",
                 position = position_dodge(0.6), width = 0.3) +
    geom_bar(aes(fill = host), alpha = 0.8, width = 0.6,
             stat = "summary", fun = "mean", position = position_dodge(0.6)) +
    geom_point(data = function(x) subset(x, group != "ref"),
               aes(group = host, color = host), size = 1, shape = 3, alpha = 0.9,
               position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1)) +
    scale_color_manual(values = point.colors) +
    scale_fill_manual(values = host.colors) +
    facet_grid(parameter~group, scales = "free", space = "free_x", switch = "y") +
    xlab("Pho4 chimera") + themes + theme(legend.position = "top")
  return(p)
}

my_plot_summary <- function(selection){
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
  selection = my_data_select(pattern = pattern, Set = c("M", "S"))
  scatter.colors = c("ScPho4" = "forestgreen", "CgPho4" = "blue3", "cyan2", "other" = "gray30")
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

my_plot_region_effect_onevar <- function(region, variable){
  # this function uses `my_calc_region_effect` output as the data
  # and makes a xy scatter plot, where x shows the value of the variable of 
  # interest with CgPho4 in the focal region, and y for the ScPho4 version
  tmp <- my_calc_region_effect(region, variable)
  p <- ggplot(tmp, aes(x = var_C, y = var_S, label = label)) +
    geom_point(size = 2.5) + 
    geom_abline(slope = 1) +
    xlab(paste0("Region ", region, " from CgPho4")) +
    ylab(paste0("Region ", region, " from ScPho4")) +
    xlim(0, NA) + ylim(0, NA) +
    ggtitle(paste0("Effect on ", variable)) +
    theme_gray(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

my_plot_region_effect_twovar_line <- function(region, highlight = "none"){
  # this function uses my_comp_region_effect to generate the data
  # and plot the difference in A_PHO2 and A_pho2 between the CgPho4 vs ScPho4
  # in the focal region
  dat <- my_comp_region_effect(region) %>% 
    pivot_longer(cols = c(dA_PHO2, dA_pho2), 
                 names_to = "host", values_to = "diff") %>% 
    mutate(host = fct_recode(host, `PHO2` = "dA_PHO2", `pho2∆` = "dA_pho2"),
           host = fct_relevel(host, "PHO2"))
  if(highlight != "none" & highlight != region){
    hl = as.numeric(highlight)
    dat <- mutate(dat, grp = str_sub(bg, hl, hl) %>% toupper(),
                  grp = fct_recode(grp, CgPho4 = "C", ScPho4 = "S"))
  }else{
    dat <- mutate(dat, grp = ifelse(NF, "n.f.", "others"))
  }
  # specify legend title
  legend.title = ""
  if(highlight != "none" & highlight != region){
    hl = as.numeric(highlight)
    dat <- mutate(dat, grp = str_sub(bg, hl, hl) %>% toupper(),
                  grp = fct_recode(grp, CgPho4 = "C", ScPho4 = "S"))
    legend.title = paste("Region", highlight, sep = " ")
  }else{
    dat <- mutate(dat, grp = ifelse(NF, "no", "yes"))
    legend.title = "Functional"
  }
  # specify arrow annotation
  arrow.x = 0.7
  arrow.y = (max(dat$diff) - min(dat$diff)) / 5 
  p <- dat %>% 
    ggplot(aes(x = host, y = diff, label = bg)) +
    geom_point(aes(color = grp), size = 2, alpha = 0.8,
               position = position_jitter(0.05)) + 
    geom_line(aes(group = bg), linewidth = 0.2, alpha = 0.8) +
    geom_segment(aes(x = arrow.x, xend = arrow.x, y = -arrow.y, yend = arrow.y),
                 arrow = arrow(length = unit(0.03, "npc"), ends = "both"), 
                 color = "gray60", lwd = 1, alpha = 0.5) +
    geom_segment(aes(x = arrow.x - 0.05, xend = arrow.x + 0.05, y = 0, yend = 0),
                 lwd = 2, color = "gray60") +
    annotate("text", x = arrow.x - 0.1, y = 3, label = "CgPho4++", 
             angle = '90', color = "gray30") +
    annotate("text", x = arrow.x + 0.1, y = -3, label = "ScPho4++", 
             angle = '270', color = "gray30") +
    scale_color_manual(legend.title, values = c("orange", "gray20")) +
    ylab("Region swap effect (Cg-Sc)") +
    theme_bw(base_size = 18) + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = rel(0.9)),
      legend.text = element_text(size = rel(0.8)),
      legend.title = element_text(size = rel(0.9)),
    )
  return(p)
}
